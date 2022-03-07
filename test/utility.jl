""" Return the number of floating point values between two values, including
one of the endpoints. """
function floatsbetween(a::T,b::T) where {T <: Float64}
  # based on https://github.com/JuliaLang/julia/commit/79920db3ac9f1e1915d836c82bc401ba0cdc322f
  # https://stackoverflow.com/questions/3587880/number-of-floats-between-two-floats
  IntT = Int64 # inttype(T)
  ia = reinterpret(IntT, a)
  ib = reinterpret(IntT, b)
  ia = ifelse(ia < zero(IntT), ia ⊻ typemax(IntT), ia)
  ib = ifelse(ib < zero(IntT), ib ⊻ typemax(IntT), ib)
  return ib-ia
end