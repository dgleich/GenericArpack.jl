macro fix_doublefloats()
  return esc( quote  
    import DoubleFloats
    import LinearAlgebra.floatmin2  
    LinearAlgebra.floatmin2(::Type{DoubleFloats.Double64}) = DoubleFloats.Double64(reinterpret(Float64, 0x2350000000000000), 0.0)    
    GenericArpack._dstqrb_maxit(::Type{DoubleFloats.Double64}) = 80
  end ) 
end 

macro fix_multifloats()
  return esc( quote
    # generate with 
    #=
    import MultiFloats
    setprecision(BigFloat, 2048)
    mftypes = [MultiFloats.Float64x2,MultiFloats.Float64x3,MultiFloats.Float64x4,MultiFloats.Float64x5,
              MultiFloats.Float64x6,MultiFloats.Float64x7,MultiFloats.Float64x8]
    println("(::Type{Int})(x::MultiFloats.MultiFloat) = Int(x._limbs[1])")
    for mft in mftypes
      println("GenericArpack._eps23(::Type{$mft}) = $(mft)(", mft(BigFloat(eps(mft))^(2/3))._limbs, ")")
      println("GenericArpack._dstqrb_maxit(::Type{$mft}) = ", 40*length(eps(mft)._limbs))
    end
    =#
    import MultiFloats
    (::Type{Int})(x::MultiFloats.MultiFloat) = Int(x._limbs[1])
    GenericArpack._eps23(::Type{MultiFloats.MultiFloat{Float64, 2}}) = MultiFloats.MultiFloat{Float64, 2}((1.3445809915232044e-21, -7.853910533714962e-38))
    GenericArpack._dstqrb_maxit(::Type{MultiFloats.MultiFloat{Float64, 2}}) = 80
    GenericArpack._eps23(::Type{MultiFloats.MultiFloat{Float64, 3}}) = MultiFloats.MultiFloat{Float64, 3}((4.930380657631343e-32, 2.3895088723650613e-49, -1.9474155911035162e-66))
    GenericArpack._dstqrb_maxit(::Type{MultiFloats.MultiFloat{Float64, 3}}) = 120
    GenericArpack._eps23(::Type{MultiFloats.MultiFloat{Float64, 4}}) = MultiFloats.MultiFloat{Float64, 4}((1.8078980427655233e-42, 5.038360808408627e-59, -4.0108760346726413e-75, 1.266024804067977e-91))
    GenericArpack._dstqrb_maxit(::Type{MultiFloats.MultiFloat{Float64, 4}}) = 160
    GenericArpack._eps23(::Type{MultiFloats.MultiFloat{Float64, 5}}) = MultiFloats.MultiFloat{Float64, 5}((6.62929611322478e-53, 1.764479965113289e-69, -7.003156249021248e-86, 1.728701256244553e-102, -8.443939604122635e-119))
    GenericArpack._dstqrb_maxit(::Type{MultiFloats.MultiFloat{Float64, 5}}) = 200
    GenericArpack._eps23(::Type{MultiFloats.MultiFloat{Float64, 6}}) = MultiFloats.MultiFloat{Float64, 6}((2.430865342914528e-63, 2.3562376651133193e-80, 4.057154919783769e-97, 1.4300476284264026e-113, -1.2991066283841956e-129, 7.263197893253085e-146))
    GenericArpack._dstqrb_maxit(::Type{MultiFloats.MultiFloat{Float64, 6}}) = 240
    GenericArpack._eps23(::Type{MultiFloats.MultiFloat{Float64, 7}}) = MultiFloats.MultiFloat{Float64, 7}((8.9136255410207e-74, -2.1331334299578177e-90, -7.205336271706034e-107, -3.334765330903578e-123, 8.258272624708863e-140, 1.1668674206531196e-156, 3.1989506894597666e-173))
    GenericArpack._dstqrb_maxit(::Type{MultiFloats.MultiFloat{Float64, 7}}) = 280
    GenericArpack._eps23(::Type{MultiFloats.MultiFloat{Float64, 8}}) = MultiFloats.MultiFloat{Float64, 8}((3.2684953330354102e-84, -9.228580946981769e-101, -4.5157222932159615e-117, 2.7992923609271783e-133, -8.868341803290226e-151, 5.0537471245946436e-167, 2.6627368329004066e-183, 1.52680335909133e-199))
    GenericArpack._dstqrb_maxit(::Type{MultiFloats.MultiFloat{Float64, 8}}) = 320
  end )
end 