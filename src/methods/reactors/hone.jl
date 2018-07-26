function hone(cur_reactor::AbstractReactor, cur_constraint::Symbol; reltol::Number=3e-3)
  cur_reactor.constraint = :beta

  prev_T = nothing
  prev_eta_CD = nothing

  max_attempts = 10
  cur_error = NaN

  for cur_index in 1:max_attempts
    println(cur_index*10)

    prev_T = cur_reactor.T_bar
    prev_eta_CD = cur_reactor.eta_CD

    new_reactor = match(cur_reactor, cur_constraint)

    cur_reactor.T_bar = new_reactor.T_bar

    tmp_reactor = deepcopy(cur_reactor)

    cur_eta_CD_list = converge(tmp_reactor; no_pts=5)

    isempty(cur_eta_CD_list) && return nothing
    @assert length(cur_eta_CD_list) == 1

    cur_bias = 0.75 / ceil(cur_index/4)
    tmp_eta_CD = cur_bias * cur_reactor.eta_CD
    tmp_eta_CD += ( 1 - cur_bias ) * cur_eta_CD_list[1]

    cur_reactor.eta_CD += cur_eta_CD_list[1]
    cur_reactor.eta_CD /= 2

    if cur_index > 1
      tmp_T_list = [prev_T, cur_reactor.T_bar]
      tmp_eta_CD_list = [prev_eta_CD, cur_reactor.eta_CD]

      cur_T_error = 1 - minimum(tmp_T_list)/maximum(tmp_T_list)
      cur_eta_CD_error = 1 - minimum(tmp_eta_CD_list)/maximum(tmp_eta_CD_list)

      cur_error = max(cur_T_error, cur_eta_CD_error)
      ( cur_error < reltol ) && break
    end
  end

  ( cur_error < reltol ) || return nothing

  honed_reactor = match(cur_reactor, cur_constraint)

  honed_reactor
end
