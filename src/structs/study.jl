mutable struct Study <: AbstractStudy
  parameter_list::Vector{AbstractFloat}

  parameter::Symbol
  sensitivity::Real
  num_points::Int

  deck::Union{Void, Symbol}

  kink_reactors::Vector{AbstractReactor}
  cost_reactors::Vector{AbstractReactor}
  wall_reactors::Vector{AbstractReactor}
end

function Study(cur_parameter; sensitivity=0.1, num_points=7, deck=nothing)
  @assert isodd(num_points)

  med_value = getfield(
    Reactor(symbols(:T_bar), deck=deck),
    cur_parameter
  )

  min_value = med_value * ( 1.0 - sensitivity )
  max_value = med_value * ( 1.0 + sensitivity )

  parameter_list = collect(linspace(min_value, max_value, num_points))

  cur_study = Study(
    parameter_list, cur_parameter,
    sensitivity, num_points, deck,
    [], [], []
  )

  for cur_value in parameter_list
    println(cur_value)
    cur_dict = Dict()
    cur_dict[:deck] = deck
    cur_dict[cur_parameter] = cur_value
    cur_dict[:constraint] = :beta

    tmp_reactor = Reactor(symbols(:T_bar), cur_dict)

    cur_kink_reactor = match(tmp_reactor, :kink)
    cur_wall_reactor = match(tmp_reactor, :wall)

    ( cur_kink_reactor == nothing ) || push!(cur_study.kink_reactors, cur_kink_reactor)
    ( cur_wall_reactor == nothing ) || push!(cur_study.wall_reactors, cur_wall_reactor)

    cur_min_T = min_T_bar
    cur_max_T = max_T_bar

    if cur_kink_reactor != nothing && cur_wall_reactor != nothing
      cur_min_T = min(cur_kink_reactor.T_bar, cur_wall_reactor.T_bar)
      cur_max_T = max(cur_kink_reactor.T_bar, cur_wall_reactor.T_bar)
    end

    cur_cost_reactor = find_min_cost_reactor(tmp_reactor, cur_min_T, cur_max_T)
    ( cur_cost_reactor == nothing ) || push!(cur_study.cost_reactors, cur_cost_reactor)
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

  @assert min_cost_reactor != nothing

  ( min_cost_reactor.T_bar == cur_T_list[1] ) && return min_cost_reactor
  ( min_cost_reactor.T_bar == cur_T_list[end] ) && return min_cost_reactor

  diff_T = ( cur_max_T - cur_min_T ) / ( num_points - 1.0 )

  beg_T = min_cost_reactor.T_bar - diff_T
  end_T = min_cost_reactor.T_bar + diff_T

  cur_func = find_min_cost(cur_reactor)

  cur_min_cost_T = Optim.minimizer(
    optimize(cur_func, beg_T, end_T; rel_tol=1e-1)
  )

  cur_reactor.T_bar = cur_min_cost_T
  work_reactor = find_low_cost_reactor(cur_reactor)

  @assert work_reactor != nothing

  min_cost_reactor = work_reactor

  min_cost_reactor
end

function find_min_cost(cur_reactor::AbstractReactor)
  cur_func = function(cur_T::Real)
    tmp_reactor = deepcopy(cur_reactor)

    tmp_reactor.T_bar = cur_T

    tmp_reactor = find_low_cost_reactor(tmp_reactor)

    ( tmp_reactor == nothing ) && return NaN

    cur_cost = tmp_reactor.cost

    cur_cost
  end

  cur_func
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
