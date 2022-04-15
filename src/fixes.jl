macro fix_doublefloats()
  return quote
    using DoubleFloats
    import LinearAlgebra.floatmin2  
    LinearAlgebra.floatmin2(::Type{Double64}) = Double64(reinterpret(Float64, 0x2350000000000000), 0.0)    
  end
end 


macro fix_multifloats()
  return quote
    GenericArpack._eps23(Float64x4) = Float64x4((1.8078980427655233e-42, 5.038360808408627e-59, -4.0108760346726413e-75, 1.2660248040679727e-91))
    GenericArpack._eps23(::Type{Float64x4}) = Float64x4((1.8078980427655233e-42, 5.038360808408627e-59, -4.0108760346726413e-75, 1.2660248040679727e-91))
    GenericArpack._dstqrb_maxit(::Type{Float64x4}) = 150
  end
end 