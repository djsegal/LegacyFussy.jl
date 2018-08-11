mutable struct Study <: AbstractStudy
  parameter_list::Vector{AbstractFloat}

  parameter::Symbol

  default::Real
  sensitivity::Real
  num_points::Int

  deck::Union{Void, Symbol}

  kink_reactors::Vector{AbstractReactor}
  cost_reactors::Vector{AbstractReactor}
  wall_reactors::Vector{AbstractReactor}
end

function Study(cur_parameter, init_dict=Dict(); sensitivity=0.2, num_points=9, verbose::Bool=true, is_parallel::Bool=true, skip_center::Bool=false, deck=nothing, cur_kwargs...)
  cur_dict = deepcopy(init_dict)

  merge!(cur_dict, Dict(cur_kwargs))

  cur_is_fast = haskey(cur_dict, :is_fast) && cur_dict[:is_fast]
  cur_is_light = haskey(cur_dict, :is_light) && cur_dict[:is_light]

  delete!(cur_dict, :is_fast)
  delete!(cur_dict, :is_light)

  @assert isodd(num_points)

  if haskey(cur_dict, :deck)
    ( deck == nothing ) &&
      ( deck = cur_dict[:deck] )
  else
    cur_dict[:deck] = deck
  end

  cur_dict[:constraint] = :beta

  cur_reactor = Reactor(symbols(:T_bar), cur_dict)

  med_value = getfield(cur_reactor, cur_parameter)

  if num_points == 1
    parameter_list = [med_value]
  else
    min_value = med_value * ( 1.0 - sensitivity )
    max_value = med_value * ( 1.0 + sensitivity )

    parameter_list = collect(linspace(min_value, max_value, num_points))
  end

  if skip_center
    deleteat!(parameter_list, 1+Int((num_points-1)/2))
    num_points -= 1
  end

  cur_study = Study(
    parameter_list, cur_parameter,
    med_value, sensitivity,
    num_points, deck,
    [], [], []
  )

  if haskey(cur_dict, :is_consistent) && cur_dict[:is_consistent]
    tmp_dict = deepcopy(cur_dict)
    delete!(tmp_dict, :is_consistent)
    other_reactor = Reactor(symbols(:T_bar), tmp_dict)

    cur_reactor.eta_CD = other_reactor.eta_CD
  end

  cur_bad_indices = []
  for (cur_index, cur_value) in enumerate(parameter_list)
    tmp_dict = deepcopy(cur_dict)
    tmp_dict[cur_parameter] = cur_value

    other_reactor = Reactor(symbols(:T_bar), tmp_dict)
    ( other_reactor.is_good ) || push!(cur_bad_indices, cur_index)
  end

  num_points -= length(cur_bad_indices)
  deleteat!(parameter_list, cur_bad_indices)

  if !cur_is_light
    cur_array = SharedArray{Float64}(num_points, 3)
    fill!(cur_array, NaN)

    cur_func = function (cur_index::Integer)
      tmp_dict = deepcopy(cur_dict)
      tmp_dict[cur_parameter] = parameter_list[cur_index]

      tmp_reactor = Reactor(symbols(:T_bar), tmp_dict)

      if tmp_reactor.is_consistent
        cur_kink_reactor = hone(tmp_reactor, :kink)
      else
        cur_kink_reactor = match(tmp_reactor, :kink)
      end

      ( cur_kink_reactor == nothing ) && return

      cur_array[cur_index,1] = cur_kink_reactor.T_bar
      cur_array[cur_index,2] = cur_kink_reactor.I_P
      cur_array[cur_index,3] = cur_kink_reactor.eta_CD
    end

    if is_parallel
      if verbose
        cur_progress = Progress(num_points)
        pmap(cur_func, cur_progress, shuffle(1:num_points))
      else
        pmap(cur_func, shuffle(1:num_points))
      end
    else
      map(cur_func, shuffle(1:num_points))
    end

    for (cur_index, cur_value) in enumerate(parameter_list)
      any(isnan, cur_array[cur_index,:]) && continue
      cur_T_bar, cur_I_P, cur_eta_CD = cur_array[cur_index,:]

      tmp_dict = deepcopy(cur_dict)
      tmp_dict[cur_parameter] = parameter_list[cur_index]

      work_reactor = Reactor(cur_T_bar, tmp_dict)
      work_reactor.I_P = cur_I_P
      work_reactor.eta_CD = cur_eta_CD

      push!(cur_study.kink_reactors, update!(work_reactor))
    end
  end

  cur_array = SharedArray{Float64}(num_points, 3)
  fill!(cur_array, NaN)

  cur_func = function (cur_index::Integer)
    tmp_dict = deepcopy(cur_dict)
    tmp_dict[cur_parameter] = parameter_list[cur_index]

    tmp_reactor = Reactor(symbols(:T_bar), tmp_dict)

    if tmp_reactor.is_consistent
      cur_wall_reactor = hone(tmp_reactor, :wall)
    else
      cur_wall_reactor = match(tmp_reactor, :wall)
    end

    ( cur_wall_reactor == nothing ) && return

    cur_array[cur_index,1] = cur_wall_reactor.T_bar
    cur_array[cur_index,2] = cur_wall_reactor.I_P
    cur_array[cur_index,3] = cur_wall_reactor.eta_CD
  end

  if is_parallel
    if verbose
      cur_progress = Progress(num_points)
      pmap(cur_func, cur_progress, shuffle(1:num_points))
    else
      pmap(cur_func, shuffle(1:num_points))
    end
  else
    map(cur_func, shuffle(1:num_points))
  end

  for (cur_index, cur_value) in enumerate(parameter_list)
    cur_T_bar, cur_I_P, cur_eta_CD = cur_array[cur_index,:]

    isnan(cur_eta_CD) && continue
    isnan(cur_T_bar) && continue
    isnan(cur_I_P) && continue

    tmp_dict = deepcopy(cur_dict)
    tmp_dict[cur_parameter] = parameter_list[cur_index]

    work_reactor = Reactor(cur_T_bar, tmp_dict)
    work_reactor.I_P = cur_I_P
    work_reactor.eta_CD = cur_eta_CD

    push!(cur_study.wall_reactors, update!(work_reactor))
  end

  if !cur_is_fast && ( cur_reactor.is_pulsed || !cur_reactor.is_consistent )
    cur_array = SharedArray{Float64}(num_points, 3)
    fill!(cur_array, NaN)

    cur_func = function (cur_index::Integer)
      tmp_dict = deepcopy(cur_dict)
      tmp_dict[cur_parameter] = parameter_list[cur_index]

      tmp_reactor = Reactor(symbols(:T_bar), tmp_dict)

      cur_kink_reactor_index = findfirst(
        tmp_kink_reactor -> getfield(tmp_kink_reactor, cur_parameter) == parameter_list[cur_index],
        cur_study.kink_reactors
      )

      cur_wall_reactor_index = findfirst(
        tmp_wall_reactor -> getfield(tmp_wall_reactor, cur_parameter) == parameter_list[cur_index],
        cur_study.wall_reactors
      )

      cur_min_T = min_T_bar
      cur_max_T = max_T_bar

      if !( iszero(cur_kink_reactor_index) || iszero(cur_wall_reactor_index) )
        cur_kink_reactor = cur_study.kink_reactors[cur_kink_reactor_index]
        cur_wall_reactor = cur_study.wall_reactors[cur_wall_reactor_index]

        cur_min_T = min(cur_kink_reactor.T_bar, cur_wall_reactor.T_bar)
        cur_max_T = max(cur_kink_reactor.T_bar, cur_wall_reactor.T_bar)
      end

      cur_cost_reactor = find_min_cost_reactor(tmp_reactor, cur_min_T, cur_max_T)
      ( cur_cost_reactor == nothing ) && return

      cur_array[cur_index,1] = cur_cost_reactor.T_bar
      cur_array[cur_index,2] = cur_cost_reactor.I_P
      cur_array[cur_index,3] = cur_cost_reactor.eta_CD
    end

    if is_parallel
      if verbose
        cur_progress = Progress(num_points)
        pmap(cur_func, cur_progress, shuffle(1:num_points))
      else
        pmap(cur_func, shuffle(1:num_points))
      end
    else
      map(cur_func, shuffle(1:num_points))
    end

    for (cur_index, cur_value) in enumerate(parameter_list)
      cur_T_bar, cur_I_P, cur_eta_CD = cur_array[cur_index,:]

      isnan(cur_eta_CD) && continue
      isnan(cur_T_bar) && continue
      isnan(cur_I_P) && continue

      tmp_dict = deepcopy(cur_dict)
      tmp_dict[cur_parameter] = parameter_list[cur_index]

      work_reactor = Reactor(cur_T_bar, tmp_dict)
      work_reactor.I_P = cur_I_P
      work_reactor.eta_CD = cur_eta_CD

      push!(cur_study.cost_reactors, update!(work_reactor))
    end
  end

  cur_study
end

function find_min_cost_reactor(cur_reactor::AbstractReactor, cur_min_T::Number, cur_max_T::Number; num_points::Int=16)
  cur_T_list = collect(linspace(cur_min_T, cur_max_T, num_points))
  out_in_reorder!(cur_T_list)

  min_cost_reactor = nothing
  is_chasing_min = [true, true]

  for (cur_index, cur_T) in enumerate(cur_T_list)
    is_chasing_min[cur_index%2+1] || @assert cur_index > 2
    is_chasing_min[cur_index%2+1] || continue

    cur_reactor.T_bar = cur_T
    work_reactor = find_low_cost_reactor(cur_reactor)

    if work_reactor == nothing
      is_chasing_min[cur_index%2+1] = false
      continue
    end

    if min_cost_reactor == nothing || work_reactor.cost < min_cost_reactor.cost
      min_cost_reactor = work_reactor
      continue
    end

    any(is_chasing_min) || break
    is_chasing_min[cur_index%2+1] = false
  end

  ( min_cost_reactor == nothing ) && return nothing

  ( min_cost_reactor.T_bar == cur_T_list[1] ) && return min_cost_reactor
  ( min_cost_reactor.T_bar == cur_T_list[end] ) && return min_cost_reactor

  diff_T = ( cur_max_T - cur_min_T ) / ( num_points - 1.0 )

  beg_T = min_cost_reactor.T_bar - diff_T
  end_T = min_cost_reactor.T_bar + diff_T

  cur_func = function(cur_T::Real)
    tmp_reactor = deepcopy(cur_reactor)
    tmp_reactor.T_bar = cur_T
    tmp_reactor = find_low_cost_reactor(tmp_reactor)

    ( tmp_reactor == nothing ) && return NaN
    tmp_reactor.cost
  end

  cur_min_cost_T = Optim.minimizer(
    optimize(cur_func, beg_T, end_T; rel_tol=1e-2)
  )

  cur_reactor.T_bar = cur_min_cost_T
  work_reactor = find_low_cost_reactor(cur_reactor)

  ( work_reactor == nothing ) && return min_cost_reactor

  min_cost_reactor = work_reactor

  min_cost_reactor
end

function find_low_cost_reactor(cur_reactor::AbstractReactor)
  cur_I_P_list = solve(cur_reactor)

  min_cost_reactor = nothing
  min_cost = Inf

  for cur_I_P in cur_I_P_list
    tmp_reactor = deepcopy(cur_reactor)
    tmp_reactor.I_P = cur_I_P
    update!(tmp_reactor)

    tmp_reactor.is_valid || continue
    ( tmp_reactor.cost < min_cost ) || continue

    min_cost_reactor = tmp_reactor
    min_cost = tmp_reactor.cost
  end

  min_cost_reactor
end
