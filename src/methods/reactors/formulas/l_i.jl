@symbol_func function l_i(cur_reactor::AbstractReactor, cur_gamma::Union{Void, Number}=nothing)
  ( cur_gamma == nothing ) &&
    ( cur_gamma = cur_reactor.gamma )

  cur_kappa = cur_reactor.kappa_95

  cur_func = function(cur_rho::Number)
    cur_numerator = ( 1.0 - cur_gamma )

    cur_denominator = 1.0

    cur_denominator += ( 1 - 3 * cur_gamma ) * cur_rho

    cur_denominator += cur_gamma * cur_rho ^ 2

    cur_value = ( cur_numerator / cur_denominator )

    cur_value ^= 4

    cur_value *= cur_rho

    cur_value
  end

  cur_l = norm_int(cur_func)

  cur_l *= 16

  cur_l *= ( 2 * cur_kappa )

  cur_l /= ( 1 + cur_kappa ^ 2 )

  cur_l
end
