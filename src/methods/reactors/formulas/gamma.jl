@symbol_func function gamma(cur_reactor::AbstractReactor, cur_rho_m::Union{Void, Number}=nothing)
  ( cur_rho_m == nothing ) &&
    ( cur_rho_m = cur_reactor.rho_m )

  cur_k_1 = 3.0 * cur_rho_m ^ 2

  cur_k_1 *= ( 3 - 4 * cur_rho_m ^ 2 + 2 * cur_rho_m ^ 4 )

  cur_k_2 = ( 3 - 6 * cur_rho_m ^ 2 + 2 * cur_rho_m ^ 4 )

  cur_k_3 = 2 - cur_rho_m ^ 2

  cur_gamma = cur_k_3

  cur_gamma /= (
    cur_k_2 + sqrt( cur_k_2 ^ 2 + cur_k_1 * cur_k_3 )
  )

  cur_gamma
end
