mutable struct Calibration <: AbstractCalibration

  default_reactor::AbstractReactor
  reactors::Vector{AbstractReactor}

  deck::Union{Void, Symbol}

  tuning_parameter::Symbol
  target_parameter::Symbol

  tuning_max_value::AbstractFloat
  target_value_list::Vector

end

function Calibration(cur_target_value_list::Vector; target_parameter::Symbol=:P_F, tuning_parameter::AbstractKey=nothing, tuning_max_value::AbstractCalculated=nothing, deck::Union{Void, Symbol}=nothing, cur_kwargs...)
  cur_dict = merge!(Dict(), Dict(cur_kwargs))
  cur_dict[:deck] = deck

  cur_reactor = Reactor(symbols(:T_bar), cur_dict)

  init_reactor = deepcopy(cur_reactor)
  _NanizeReactor!(init_reactor)
  init_reactor.T_bar = NaN

  ( tuning_parameter == nothing ) &&
    ( tuning_parameter = cur_reactor.is_pulsed ? :B_CS : :H )

  ( tuning_max_value == nothing ) &&
    ( tuning_max_value = getfield(cur_reactor, tuning_parameter) * 3 )

  cur_calibration = Calibration(
    init_reactor, [], deck,
    tuning_parameter, target_parameter,
    tuning_max_value, cur_target_value_list
  )

  num_points = length(cur_target_value_list)

  cur_array = SharedArray{Float64}(num_points, 4)
  fill!(cur_array, NaN)

  cur_func = function (cur_index)
    cur_target_value = cur_target_value_list[cur_index]

    cur_prototype = find_prototype(
      cur_reactor, tuning_parameter, target_parameter,
      tuning_max_value, cur_target_value
    )

    ( cur_prototype == nothing ) && return

    cur_vector = [
      cur_prototype.T_bar,
      cur_prototype.I_P,
      cur_prototype.eta_CD,
      getfield(cur_prototype, tuning_parameter)
    ]

    cur_array[cur_index,:] = cur_vector

    return
  end

  cur_progress = Progress(num_points)
  pmap(cur_func, cur_progress, shuffle(1:num_points))

  for cur_index in 1:num_points
    cur_vector = cur_array[cur_index,:]

    all(!isnan, cur_vector) || continue
    cur_T_bar, cur_I_P, cur_eta_CD, cur_tuned_value = cur_vector

    tmp_reactor = deepcopy(cur_reactor)

    tmp_reactor.T_bar = cur_T_bar
    tmp_reactor.I_P = cur_I_P
    tmp_reactor.eta_CD = cur_eta_CD

    setfield!(tmp_reactor, tuning_parameter, cur_tuned_value)
    push!(cur_calibration.reactors, update!(tmp_reactor))
  end

  cur_calibration

end


function find_prototype(cur_reactor::AbstractReactor, tuning_parameter::Symbol, target_parameter::Symbol, tuning_max_value::AbstractFloat, cur_target_value::AbstractFloat)
  tmp_func = function(tmp_value)
    cur_value = real(tmp_value)
    isnan(cur_value) && return NaN

    tmp_reactor = deepcopy(cur_reactor)
    setfield!(tmp_reactor, tuning_parameter, cur_value)

    if tmp_reactor.is_consistent
      tmp_reactor = hone(tmp_reactor, :wall)
    else
      tmp_reactor = match(tmp_reactor, :wall)
    end

    ( tmp_reactor == nothing ) && return NaN
    tmp_reactor.is_valid || return NaN

    cur_error = cur_target_value - getfield(tmp_reactor, target_parameter)
    cur_error
  end

  cur_root_max = getfield(cur_reactor, tuning_parameter) * 3.0
  cur_roots = find_roots(tmp_func, 0.0, tuning_max_value, abstol=1e-4, no_pts = 11)

  isempty(cur_roots) && return nothing
  @assert length(cur_roots) == 1
  cur_root = cur_roots[1]

  cur_prototype = deepcopy(cur_reactor)
  setfield!(cur_prototype, tuning_parameter, cur_root)

  if cur_prototype.is_consistent
    cur_prototype = hone(cur_prototype, :wall)
  else
    cur_prototype = match(cur_prototype, :wall)
  end

  cur_error = cur_target_value - cur_prototype.P_F
  @assert cur_error < 0.01

  cur_prototype
end
