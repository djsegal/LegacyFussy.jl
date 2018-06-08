@symbol_func function perimeter(cur_reactor::AbstractReactor)
  cur_perimeter = 2.0

  cur_perimeter *= cur_reactor.pi

  cur_perimeter *= a(cur_reactor)

  cur_perimeter *= sqrt(
    _perimeter_g(cur_reactor)
  )

  cur_perimeter
end
