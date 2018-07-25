mutable struct Study <: AbstractStudy
  T_bar_list::Vector{AbstractFloat}

  parameter::Symbol
  sensitivity::Real
  num_points::Int

  deck::Union{Void, Symbol}

  kink_reactors::Vector{AbstractReactor}
  cost_reactors::Vector{AbstractReactor}
  wall_reactors::Vector{AbstractReactor}
end

function Study(cur_T_bar_list, cur_parameter; sensitivity=0.1, num_points=7, deck=nothing)
  cur_study = Study(
    collect(cur_T_bar_list), cur_parameter,
    sensitivity, num_points, deck,
    [], [], []
  )

  med_value = getfield(
    Reactor(first(cur_T_bar_list), deck=deck),
    cur_parameter
  )

  min_value = med_value * ( 1.0 - sensitivity )
  max_value = med_value * ( 1.0 + sensitivity )

  cur_values = collect(linspace(min_value, max_value, num_points))

  cur_scans = []

  @showprogress 1 "Progress: " for cur_value in cur_values
    cur_dict = Dict()
    cur_dict[:deck] = deck
    cur_dict[cur_parameter] = cur_value
    cur_dict[:ignored_limits] = [:kink, :wall, :heat]
    cur_dict[:max_I_P] = 40
    cur_dict[:no_pts_I_P] = 201
    cur_dict[:verbose] = false

    push!(cur_scans, Scan(cur_T_bar_list, :beta; cur_dict...))
  end

  cur_array = SharedArray{Float64}(num_points, 3)
  fill!(cur_array, NaN)

  cur_func = function (cur_index::Integer)
    cur_array[cur_index, :] = study_scan(cur_scans[cur_index], cur_parameter)
  end

  cur_progress = Progress(num_points)
  pmap(cur_func, cur_progress, shuffle(1:num_points))

  cur_scan_reactors = [
    cur_study.kink_reactors,
    cur_study.cost_reactors,
    cur_study.wall_reactors
  ]

  for (cur_index, cur_scan) in enumerate(cur_scans)
    for (cur_sub_index, cur_T) in enumerate(cur_array[cur_index,:])
      isnan(cur_T) && continue

      tmp_reactor = deepcopy(cur_scan.beta_reactors[1])
      tmp_reactor.T_bar = cur_T

      cur_reactor_list = _get_raw_reactor_list(tmp_reactor, cur_parameter, [:heat])
      filter!(cur_reactor -> cur_reactor.is_valid, cur_reactor_list)

      @assert 1 <= cur_sub_index <= 3

      if cur_sub_index == 2
        cur_min_index = indmin(map(cur_reactor -> cur_reactor.cost, cur_reactor_list))
      elseif cur_sub_index == 1
        cur_min_index = indmin(map(cur_reactor -> abs( 1 - cur_reactor.norm_q_95 ), cur_reactor_list))
      elseif cur_sub_index == 3
        cur_min_index = indmin(map(cur_reactor -> abs( 1 - cur_reactor.norm_P_W ), cur_reactor_list))
      end

      push!(cur_scan_reactors[cur_sub_index], cur_reactor_list[cur_min_index])
    end
  end

  cur_study
end

function study_scan(cur_scan::AbstractScan, cur_parameter::Symbol)
  (kink_T, wall_T, min_T_list, max_T_list) = find_scan_borders(cur_scan, cur_parameter)

  ( isempty(min_T_list) || isempty(max_T_list) ) && return [NaN, NaN, NaN]

  cost_T = find_scan_min(cur_scan, cur_parameter, min_T_list, max_T_list)

  [kink_T, cost_T, wall_T]
end

function _process_scan_reactors!(cur_scan::AbstractScan, cur_parameter::Symbol)
  cur_reactor_list = filter(cur_reactor -> cur_reactor.is_valid, cur_scan.beta_reactors)
  isempty(cur_reactor_list) && return []

  sort!(cur_reactor_list, by = cur_reactor -> cur_reactor.T_bar)

  reactor_count = Int( 1 + ceil(length(cur_scan.beta_reactors)/4) )

  beg_reactors = []
  end_reactors = []

  tmp_reactor = cur_reactor_list[1]

  cur_dict = Dict()
  cur_dict[:deck] = tmp_reactor.deck
  cur_dict[cur_parameter] = getfield(tmp_reactor, cur_parameter)
  cur_dict[:ignored_limits] = [:kink, :wall, :heat]
  cur_dict[:max_I_P] = 40
  cur_dict[:no_pts_I_P] = 201
  cur_dict[:is_parallel] = false

  min_T = first(cur_reactor_list).T_bar
  max_T = last(cur_reactor_list).T_bar

  if min_T != minimum(cur_scan.T_bar_list)
    prev_T = sort(cur_scan.T_bar_list)[findfirst(cur_T -> cur_T == min_T, cur_scan.T_bar_list) - 1]
    tmp_T_bar_list = collect(linspace(prev_T, min_T, reactor_count))[2:end]

    beg_reactors = Scan(tmp_T_bar_list, :beta; cur_dict...).beta_reactors
    filter!(cur_reactor -> cur_reactor.is_valid, beg_reactors)
    sort!(beg_reactors, by = cur_reactor -> cur_reactor.T_bar)

    this_reactor = beg_reactors[end]
    that_reactor = cur_reactor_list[1]

    @assert this_reactor.T_bar == min_T
    @assert that_reactor.T_bar == min_T

    cur_error = abs( this_reactor.R_0 - that_reactor.R_0 )
    for work_reactor in cur_reactor_list[2:end]
      ( work_reactor.T_bar == min_T ) || break

      tmp_error = abs( this_reactor.R_0 - work_reactor.R_0 )
      ( tmp_error < cur_error ) || continue

      cur_error = tmp_error
      that_reactor = work_reactor
    end

    @assert isapprox(cur_error, 0.0, atol=1e-4)

    cur_shift = ( that_reactor.branch_id - this_reactor.branch_id )

    if !iszero(cur_shift)
      for work_reactor in beg_reactors
        work_reactor.branch_id += cur_shift
      end
    end

    pop!(beg_reactors)
  end

  if max_T != maximum(cur_scan.T_bar_list)
    next_T = sort(cur_scan.T_bar_list)[findfirst(cur_T -> cur_T == max_T, cur_scan.T_bar_list) + 1]
    tmp_T_bar_list = collect(linspace(max_T, next_T, reactor_count))[1:end-1]

    end_reactors = Scan(tmp_T_bar_list, :beta; cur_dict...).beta_reactors
    woof = deepcopy(end_reactors)
    filter!(cur_reactor -> cur_reactor.is_valid, end_reactors)
    sort!(end_reactors, by = cur_reactor -> cur_reactor.T_bar)

    this_reactor = end_reactors[1]
    that_reactor = cur_reactor_list[end]

    if this_reactor.T_bar != max_T || that_reactor.T_bar != max_T
      println(111, " - ", this_reactor)
      println(222, " - ", end_reactors)
      println(202, " - ", woof)
      println(333, " - ", cur_reactor_list)
      println(444, " - ", max_T)
    end

    @assert this_reactor.T_bar == max_T
    @assert that_reactor.T_bar == max_T

    cur_error = abs( this_reactor.R_0 - that_reactor.R_0 )
    for work_reactor in reverse(cur_reactor_list)[2:end]
      ( work_reactor.T_bar == max_T ) || break

      tmp_error = abs( this_reactor.R_0 - work_reactor.R_0 )
      ( tmp_error < cur_error ) || continue

      cur_error = tmp_error
      that_reactor = work_reactor
    end

    @assert isapprox(cur_error, 0.0, atol=1e-4)

    cur_shift = ( that_reactor.branch_id - this_reactor.branch_id )

    if !iszero(cur_shift)
      for work_reactor in end_reactors
        work_reactor.branch_id += cur_shift
      end
    end

    shift!(end_reactors)
  end

  prepend!(cur_reactor_list, beg_reactors)
  append!(cur_reactor_list, end_reactors)

  cur_branch_list = map(cur_reactor -> cur_reactor.branch_id, cur_reactor_list)
  uniq_branch_list = unique(cur_branch_list)

  min_branch_id = minimum(cur_branch_list)
  max_branch_id = maximum(cur_branch_list)

  if min_branch_id != 1
    cur_shift = ( 1 - min_branch_id )
    max_branch_id += cur_shift

    for work_reactor in cur_reactor_list
      work_reactor.branch_id += cur_shift
    end
  end

  @assert max_branch_id == length(uniq_branch_list)
  @assert max_branch_id < 3

  cur_reactor_list
end

function find_scan_borders(cur_scan::Scan, cur_parameter::Symbol)
  cur_reactor_list = _process_scan_reactors!(cur_scan, cur_parameter)

  valid_reactors = filter(
    cur_reactor -> ( cur_reactor.norm_q_95 <= 1 && cur_reactor.norm_P_W <= 1 ),
    cur_reactor_list
  )

  iszero(length(valid_reactors)) && return [NaN, NaN, [], []]

  min_T_list = []
  max_T_list = []

  max_branch_id = maximum(map(cur_reactor -> cur_reactor.branch_id, cur_reactor_list))

  for cur_branch_id in 1:max_branch_id
    work_T_bar_list = map(
      cur_reactor -> cur_reactor.T_bar,
      filter(work_reactor -> ( work_reactor.branch_id == cur_branch_id ), valid_reactors)
    )

    if isempty(work_T_bar_list)
      push!(min_T_list, NaN)
      push!(max_T_list, NaN)
    else
      push!(min_T_list, minimum(work_T_bar_list))
      push!(max_T_list, maximum(work_T_bar_list))
    end
  end

  bad_kink_reactors = filter(
    cur_reactor -> ( cur_reactor.norm_q_95 > 1 && cur_reactor.norm_P_W <= 1 ),
    cur_reactor_list
  )

  bad_wall_reactors = filter(
    cur_reactor -> ( cur_reactor.norm_q_95 <= 1 && cur_reactor.norm_P_W > 1 ),
    cur_reactor_list
  )

  valid_branch_list = map(cur_reactor -> cur_reactor.branch_id, valid_reactors)

  if !isempty(bad_kink_reactors)
    bad_kink_branch_list = map(cur_reactor -> cur_reactor.branch_id, bad_kink_reactors)
    filter!(tmp_branch_id -> in(tmp_branch_id, valid_branch_list), bad_kink_branch_list)

    if isempty(bad_kink_branch_list)
      empty!(bad_kink_reactors)
    else
      ( length(unique(bad_kink_branch_list)) == 1 ) || println("~~~" , bad_kink_reactors)
      @assert length(unique(bad_kink_branch_list)) == 1
      bad_kink_branch_id = bad_kink_branch_list[1]
      filter!(cur_reactor -> cur_reactor.branch_id == bad_kink_branch_id, bad_kink_reactors)

      bad_kink_T_bar_list = map(cur_reactor -> cur_reactor.T_bar, bad_kink_reactors)
      min_bad_kink_T = minimum(bad_kink_T_bar_list)
      max_bad_kink_T = maximum(bad_kink_T_bar_list)
    end
  end

  if !isempty(bad_wall_reactors)
    bad_wall_branch_list = map(cur_reactor -> cur_reactor.branch_id, bad_wall_reactors)
    filter!(tmp_branch_id -> in(tmp_branch_id, valid_branch_list), bad_wall_branch_list)

    if isempty(bad_wall_branch_list)
      empty!(bad_wall_reactors)
    else
      ( length(unique(bad_wall_branch_list)) == 1 ) || println("___" , bad_wall_reactors)
      @assert length(unique(bad_wall_branch_list)) == 1
      bad_wall_branch_id = bad_wall_branch_list[1]
      filter!(cur_reactor -> cur_reactor.branch_id == bad_wall_branch_id, bad_wall_reactors)

      bad_wall_T_bar_list = map(cur_reactor -> cur_reactor.T_bar, bad_wall_reactors)
      min_bad_wall_T = minimum(bad_wall_T_bar_list)
      max_bad_wall_T = maximum(bad_wall_T_bar_list)
    end
  end

  kink_T = NaN
  wall_T = NaN

  if !isempty(bad_kink_reactors)
    kink_T_func = find_beta_kink_intersection(valid_reactors[1], cur_parameter, bad_kink_branch_id)

    if !isnan(min_T_list[bad_kink_branch_id])
      if ( max_bad_kink_T < min_T_list[bad_kink_branch_id] )
        kink_T_guesses = (max_bad_kink_T, min_T_list[bad_kink_branch_id])
        kink_T = find_zero(kink_T_func, kink_T_guesses, FalsePosition())
        min_T_list[bad_kink_branch_id] = kink_T
      else
        kink_T_guesses = (max_T_list[bad_kink_branch_id], min_bad_kink_T)
        kink_T = find_zero(kink_T_func, kink_T_guesses, FalsePosition())
        max_T_list[bad_kink_branch_id] = kink_T
      end
    end
  end

  if !isempty(bad_wall_reactors)
    wall_T_func = find_beta_wall_intersection(valid_reactors[1], cur_parameter, bad_wall_branch_id)

    if !isnan(min_T_list[bad_wall_branch_id])
      if ( max_bad_wall_T < min_T_list[bad_wall_branch_id] )
        wall_T_guesses = (max_bad_wall_T, min_T_list[bad_wall_branch_id])
        wall_T = find_zero(wall_T_func, wall_T_guesses, FalsePosition())
        min_T_list[bad_wall_branch_id] = wall_T
      else
        wall_T_guesses = (max_T_list[bad_wall_branch_id], min_bad_wall_T)
        wall_T = find_zero(wall_T_func, wall_T_guesses, FalsePosition())
        max_T_list[bad_wall_branch_id] = wall_T
      end
    end
  end

  return (kink_T, wall_T, min_T_list, max_T_list)
end

function find_scan_min(cur_scan::Scan, cur_parameter::Symbol, min_T_list::Vector, max_T_list::Vector)
  sample_reactor = first(cur_scan.beta_reactors)

  work_T_list = []

  for (cur_branch_id, (cur_min_T, cur_max_T)) in enumerate(zip(min_T_list, max_T_list))
    ( isnan(cur_min_T) || isnan(cur_max_T) ) && continue

    cur_func = find_min_cost(sample_reactor, cur_parameter, cur_branch_id)

    push!(work_T_list, Optim.minimizer(optimize(cur_func, cur_min_T, cur_max_T, rel_tol=3e-3)))
  end

  cur_cost = Inf
  cost_T = NaN

  append!(work_T_list, min_T_list)
  append!(work_T_list, max_T_list)

  work_T_list = unique(work_T_list)
  filter!(!isnan, work_T_list)

  for cur_T in unique([ work_T_list..., min_T_list... , max_T_list... ])
    tmp_reactor = deepcopy(sample_reactor)
    tmp_reactor.T_bar = cur_T
    cur_reactor_list = _get_raw_reactor_list(tmp_reactor, cur_parameter, [:heat])

    filter!(cur_reactor -> cur_reactor.is_valid, cur_reactor_list)
    isempty(cur_reactor_list) && continue

    tmp_cost = minimum(map(work_reactor -> work_reactor.cost, cur_reactor_list))
    ( tmp_cost < cur_cost ) || continue

    cur_cost = tmp_cost
    cost_T = cur_T
  end

  cost_T
end

function find_min_cost(cur_reactor::AbstractReactor, cur_parameter::Symbol, cur_branch_id::Int)
  cur_func = function(cur_T::Real)
    tmp_reactor = deepcopy(cur_reactor)

    tmp_reactor.T_bar = cur_T

    cur_reactor_list = _get_raw_reactor_list(tmp_reactor, cur_parameter, [:heat])

    filter!(
      cur_reactor -> cur_reactor.is_valid,
      cur_reactor_list
    )

    work_reactor = cur_reactor_list[min(length(cur_reactor_list),abs(cur_branch_id))]

    cur_cost = work_reactor.cost

    cur_cost
  end

  cur_func
end

function find_beta_wall_intersection(cur_reactor::AbstractReactor, cur_parameter::Symbol, cur_branch_id::Int)
  cur_func = function(cur_T::Real)
    tmp_reactor = deepcopy(cur_reactor)

    tmp_reactor.T_bar = cur_T

    cur_reactor_list = _get_raw_reactor_list(tmp_reactor, cur_parameter, [:heat, :wall])

    filter!(
      cur_reactor -> cur_reactor.is_valid,
      cur_reactor_list
    )

    work_reactor = cur_reactor_list[min(length(cur_reactor_list),abs(cur_branch_id))]

    cur_error = ( 1.0 - work_reactor.norm_P_W )

    cur_error
  end

  cur_func
end

function find_beta_kink_intersection(cur_reactor::AbstractReactor, cur_parameter::Symbol, cur_branch_id::Int)
  cur_func = function(cur_T::Real)
    tmp_reactor = deepcopy(cur_reactor)

    tmp_reactor.T_bar = cur_T

    cur_reactor_list = _get_raw_reactor_list(tmp_reactor, cur_parameter, [:heat, :kink])

    filter!(
      cur_reactor -> cur_reactor.is_valid,
      cur_reactor_list
    )

    work_reactor = cur_reactor_list[min(length(cur_reactor_list),abs(cur_branch_id))]

    cur_error = ( 1.0 - work_reactor.norm_q_95 )

    cur_error
  end

  cur_func
end

function _get_raw_reactor_list(cur_reactor::AbstractReactor, cur_parameter::Symbol, ignored_limits::Vector)
  tmp_reactor = Reactor(
    cur_reactor.T_bar,
    deck = cur_reactor.deck,
    constraint = :beta,
    ignored_limits = ignored_limits
  )

  setfield!(
    tmp_reactor, cur_parameter,
    getfield(cur_reactor, cur_parameter)
  )

  cur_I_P_list = solve(tmp_reactor)

  cur_reactor_list = []

  for cur_I_P in cur_I_P_list
    work_reactor = deepcopy(tmp_reactor)
    work_reactor.I_P = cur_I_P

    push!(cur_reactor_list, update!(work_reactor))
  end

  cur_reactor_list
end
