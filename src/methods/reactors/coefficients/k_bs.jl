@symbol_func function K_BS(cur_reactor::AbstractReactor)
  cur_nu_n = cur_reactor.nu_n

  cur_nu_T = cur_reactor.nu_T

  cur_K = 0.6099

  cur_K *= ( 1 + cur_reactor.kappa_95 ^ 2 )

  cur_K /= 2

  cur_K *= sqrt( cur_reactor.epsilon ) ^ 5

  cur_K *= ( 1 + cur_nu_n )

  cur_K *= ( 1 + cur_nu_T )

  cur_K *= ( cur_nu_n + 0.054 * cur_nu_T )

  cur_K *= _C_BS(cur_reactor)

  cur_K
end

@symbol_func function _C_BS(cur_reactor::AbstractReactor)
  cur_nu_n = cur_reactor.nu_n

  cur_nu_T = cur_reactor.nu_T

  cur_gamma = cur_reactor.gamma

  cur_func = function (cur_rho)
    cur_value = ( 1.0 - cur_rho )

    cur_value ^= ( cur_nu_n + cur_nu_T - 1 )

    cur_value *= cur_rho ^ (1/4)

    cur_value *= (
      1.0 +
      ( 1 - 3 * cur_gamma ) * cur_rho +
      ( cur_gamma ) * cur_rho ^ 2
    ) ^ 2

    cur_value
  end

  cur_reactor.is_symbolic &&
    ( cur_func = cur_func(rho_sym) )

  cur_C = norm_int(cur_func)

  cur_C /= ( 1 - cur_gamma ) ^ 2

  cur_C
end
