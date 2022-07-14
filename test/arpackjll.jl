import Arpack_jll, LinearAlgebra, SHA

# show the library if we are debugging...
@debug Arpack_jll.libarpack, bytes2hex(open(SHA.sha256, Arpack_jll.libarpack))

function _reset_libarpack_dgetv0_iseed()
  # we need hardcoded offsets because we can't get these from cglobal :( )
  sha = open(SHA.sha256, Arpack_jll.libarpack_path)
  if bytes2hex(sha) == "5c1e951fad68bd7b180b83f2d821324efbdda7e00f5926698fad720249d6ac3f"
    # iseedoffsets = (0x0000000000059fc0, ) # this is for just dgetv0 
    # this is for all of them! 
    iseedoffsets = (0x0000000000059fc0, 0x0000000000059fe0, 0x0000000000059fa0, 0x000000000005a000)
    dgetv0offset = 0x000000000001f940
  elseif bytes2hex(sha) == "84b5fde039119bd1f80489900b57f7fcf272af081bb37fae7d3a90e6ca673569"
    dgetv0offset = 0x000000000001c570
    iseedoffsets = (0x00000000000538c0, 0x0000000000053d20, 0x0000000000053680, 0x0000000000054120)
  elseif bytes2hex(sha) == "88d0208af594702e43049e01e813d6ca5e15cc4974a9db782fbf74fe7cb93613"
    dgetv0offset = 0x000000000001cb60
    iseedoffsets = (0x00000000000538c0, 0x0000000000053d20, 0x0000000000053680, 0x0000000000054120)
  elseif bytes2hex(sha) == "db784b23887f4de4e9ea184e252c8efa356c4a96b5f658298627fb6dd237a3f3"
    dgetv0offset = 0x000000000001e4f0
    iseedoffsets = (0x0000000000058fc0, 0x0000000000058fe0, 0x0000000000058fa0, 0x0000000000059000)
  elseif bytes2hex(sha) == "9dc878a0bdcce61ddbd725908e592ad4d4e104922d6583158599cac6b0b5e029"
    dgetv0offset = 0x000000000001b4c0
    iseedoffsets = (0x00000000002538c0, 0x0000000000253d20, 0x0000000000253680, 0x0000000000254120)
  elseif bytes2hex(sha) == "ac654701c6be6af0fc57c05f5d3c60523c66a5c76e7940726b0316504bf7e9a8"
    dgetv0offset = 0x000000000001b800
    iseedoffsets = (0x00000000002538c0, 0x0000000000253d20, 0x0000000000253680, 0x0000000000254120)
  elseif bytes2hex(sha) == "fdd499ed4c0826f147ca4d57f4b0fad9edebc8aae6a5b038c1a57f8647b515b7"
    dgetv0offset = 0x000000000001e7b0
    iseedoffsets = (0x000000000005afc0, 0x000000000005afe0, 0x000000000005afa0, 0x000000000005b000)  
  elseif bytes2hex(sha) == "8b5b5e500dc839928a5ce4f124a84afa8331d99b22775211ec098c5bf4280c6f"
    # this is from windows 
    # Quick study of windows format to get the right offsets here...
    # 
    #dgetv0offset = 0x19400
    # objdump -t -p ~/Downloads/Arpack.v3.5.1.x86_64-w64-mingw32-libgfortran5/bin/libarpack.dll | grep dgetv0    
    #       15  0x19400  dgetv0_
    #  AUX dgetv0.f
    #  [581](sec  1)(fl 0x00)(ty  20)(scl   2) (nx 1) 0x00018400 dgetv0_
    # we use the second one, because that's the more general symbol offset
    # [326](sec  6)(fl 0x00)(ty   0)(scl   3) (nx 0) 0x00000060 iseed.3897
    # [584](sec  6)(fl 0x00)(ty   0)(scl   3) (nx 0) 0x000002c0 iseed.3896
    # [1147](sec  6)(fl 0x00)(ty   0)(scl   3) (nx 0) 0x00000780 iseed.3896
    # [1710](sec  6)(fl 0x00)(ty   0)(scl   3) (nx 0) 0x00000be0 iseed.3897
    # above, we have relative offsets to the second base... ugh, annoying! 

    # these are almost correct... 
    #dgetv0offset = 0x00019400
    #iseedoffsets = (0x00000060, 0x000002c0, 0x00000780, 0x00000be0)

    # look at the disassembly
    # dgleich@circulant GenericArpack % objdump -td -p ~/Downloads/Arpack.v3.5.1.x86_64-w64-mingw32-libgfortran5/bin/libarpack.dll | grep iseed 
    # 64510a86: 0f 29 05 d3 45 04 00         	movaps	%xmm0, 280019(%rip)     # 0x64555060 <iseed.3897>
    # 64510a95: 0f 29 05 d4 45 04 00         	movaps	%xmm0, 280020(%rip)     # 0x64555070 <iseed.3897+0x10>
    # 64511078: 48 8d 15 e1 3f 04 00         	leaq	278497(%rip), %rdx      # 0x64555060 <iseed.3897>
    # 64519453: 0f 29 05 66 be 03 00         	movaps	%xmm0, 245350(%rip)     # 0x645552c0 <iseed.3896>
    # 64519462: 0f 29 05 67 be 03 00         	movaps	%xmm0, 245351(%rip)     # 0x645552d0 <iseed.3896+0x10>
    # 645199e0: 48 8d 15 d9 b8 03 00         	leaq	243929(%rip), %rdx      # 0x645552c0 <iseed.3896>
    # 6452a193: 0f 29 05 e6 b5 02 00         	movaps	%xmm0, 177638(%rip)     # 0x64555780 <iseed.3896>
    # 6452a1a2: 0f 29 05 e7 b5 02 00         	movaps	%xmm0, 177639(%rip)     # 0x64555790 <iseed.3896+0x10>
    # 6452a710: 48 8d 15 69 b0 02 00         	leaq	176233(%rip), %rdx      # 0x64555780 <iseed.3896>
    # 6453b646: 0f 29 05 93 a5 01 00         	movaps	%xmm0, 107923(%rip)     # 0x64555be0 <iseed.3897>
    # 6453b655: 0f 29 05 94 a5 01 00         	movaps	%xmm0, 107924(%rip)     # 0x64555bf0 <iseed.3897+0x10>
    # 6453bc50: 48 8d 15 89 9f 01 00         	leaq	106377(%rip), %rdx      # 0x64555be0 <iseed.3897>  
    
    iseedoffsets = (0x64555060, 0x645552c0, 0x64555be0, 0x64555780)
    dgetv0offset = 0x64519400

    # note that these are consistent 
    # dgetv0offset = 0x64519400
    #                     ^^^^^
    # iseedoffsets = (0x64555060, 0x645552c0, 0x64555be0, 0x64555780)
    #                      ^^^^^         ^^^         ^^^         ^^^
    
    # so somewehre, we must get: section 1 has offset 0x645000000 and 
    # this is from __ImageBase ? 
    # [3677](sec -1)(fl 0x00)(ty   0)(scl   2) (nx 0) 0x64500000 __ImageBase
  elseif bytes2hex(sha) == "8c3e31745da5072b279e9c4b6c5873b2a92d0c431e9c98fe7b7e6ae682cbf0ae"
    # Another windows case...

    
    iseedoffsets = (0x64554be0, 0x64554780, 0x645542c0, 0x64554060)
    dgetv0offset = 0x64518da0
  else
    # try and get them automatically...
    if Sys.isapple() || Sys.islinux()
      run(pipeline(`nm $(Arpack_jll.libarpack_path)`,`grep dgetv0`))
      run(pipeline(`nm $(Arpack_jll.libarpack_path)`,`grep iseed`))
      @show bytes2hex(open(SHA.sha256, Arpack_jll.libarpack_path))
      iseed_offsets_str = readchomp(pipeline(`nm $(Arpack_jll.libarpack_path)`,`grep iseed`))
      iseed_offsets = iseed_offsets_str |> 
                x->split(x,"\n") |>
                x->map(y->y[1:16],x) |>
                x->parse.(UInt64, x;base=16)
      dgetv0_offset = readchomp(pipeline(`nm $(Arpack_jll.libarpack_path)`,`grep dgetv0`)) |> x->parse(UInt64, x[1:16];base=16)

      errmsg = """Unknown libarpack sha $(bytes2hex(sha)) for $(Arpack_jll.libarpack_path)

      got iseed_offsets= $(iseed_offsets)

      got dgetv0_offset = $(repr(dgetv0_offset))

      Please post a new issue on the github page for GenericArpack."""
      @error(errmsg)
    elseif Sys.iswindows()
      run(pipeline(`objdump -td -p $(Arpack_jll.libarpack_path)`,`grep dgetv0`))
      run(pipeline(`objdump -td -p $(Arpack_jll.libarpack_path)`,`grep iseed`))
      run(pipeline(`objdump -td -p $(Arpack_jll.libarpack_path)`,`grep __ImageBase`))
      errmsg = """Unknown libarpack sha $(bytes2hex(sha)) for $(Arpack_jll.libarpack_path)

      Please post a new issue on the github page for GenericArpack."""
      @error(errmsg)
    else

      errmsg = """Unknown libarpack sha $(bytes2hex(sha)) for $(Arpack_jll.libarpack_path)

      Please post a new issue on the github page for GenericArpack."""
      @error(errmsg)
    end 
  end
  libar = Base.Libc.dlopen(Arpack_jll.libarpack_path)
  dgetv0_real_offset = Base.Libc.dlsym(libar, "dgetv0_")
  base_offset = dgetv0_real_offset-dgetv0offset # this comes from the command above

  for iseedoffset in iseedoffsets 
    piseedoffset = Ptr{LinearAlgebra.BlasInt}(base_offset + iseedoffset)
    #previseed = unsafe_load.(piseedoffset, (1,2,3,4))


    # store the values 1,3,5,7 at indices 1,2,3,4... which resets to the initial Arpack config.
    unsafe_store!.(piseedoffset, (1,3,5,7), (1,2,3,4))
    #newiseed = unsafe_load.(piseedoffset, (1,2,3,4))
  end 

  #return previseed, newiseed
  return true
end 

function arpack_set_debug_high()
  # [Documentation and structure of debug block here.](https://github.com/opencollab/arpack-ng/blob/master/SRC/debug.h)
  arpack_debug = cglobal((:debug_, Arpack_jll.libarpack), Int64)
  unsafe_store!(arpack_debug, 6, 1) # logfile set logfil to 6, the default stdout
  unsafe_store!(arpack_debug, -20, 2) # ndigit - use 14 digits of precision (ndigit)
  unsafe_store!.(arpack_debug, 4, 3:24) # turn on most debugging
end

function arpack_set_debug_low()
  # [Documentation and structure of debug block here.](https://github.com/opencollab/arpack-ng/blob/master/SRC/debug.h)
  arpack_debug = cglobal((:debug_, Arpack_jll.libarpack), Int64)
  unsafe_store!(arpack_debug, 6, 1) # logfile set logfil to 6, the default stdout
  unsafe_store!(arpack_debug, -6, 2) # ndigit - use 14 digits of precision (ndigit)
  unsafe_store!.(arpack_debug, 1, 3:24) # turn on most debugging
end

function arpack_set_debug_off()
  # [Documentation and structure of debug block here.](https://github.com/opencollab/arpack-ng/blob/master/SRC/debug.h)
  arpack_debug = cglobal((:debug_, Arpack_jll.libarpack), Int64)
  unsafe_store!(arpack_debug, 6, 1) # logfile set logfil to 6, the default stdout
  unsafe_store!(arpack_debug, -6, 2) # ndigit - use 14 digits of precision (ndigit)
  unsafe_store!.(arpack_debug, 0, 3:24) # turn off most debugging
end

function arpack_aupd_output()
  # [Documentation and structure of debug block here.](https://github.com/opencollab/arpack-ng/blob/master/SRC/debug.h)
  arpack_debug = cglobal((:debug_, Arpack_jll.libarpack), Int64)
  unsafe_store!(arpack_debug, 6, 1) # logfile set logfil to 6, the default stdout
  unsafe_store!(arpack_debug, -6, 2) # ndigit - use 14 digits of precision (ndigit)
  unsafe_store!.(arpack_debug, 1, [4 11 18]) # turn on aupd output 
end



function arpack_dsconv(n::Int, ritz::Vector{Float64}, bounds::Vector{Float64},
                      tol::Float64)
  nconv = Ref{LinearAlgebra.BlasInt}(-1)
  ccall((:dsconv_, Arpack_jll.libarpack), Cvoid,
    (Ref{LinearAlgebra.BlasInt},
     Ptr{Float64},
     Ptr{Float64},
     Ref{Float64},
     Ref{LinearAlgebra.BlasInt}),
    n, ritz, bounds, tol, nconv)
  return nconv[]
end

function arpack_dsortr(
  which::Symbol, # Input
  apply::Bool, # Input
  n::Int, # Input
  x1::Vector{Float64}, # Input/Output
  x2::Vector{Float64}, # Input/Output
  )
  whichstr = string(which)
  ccall((:dsortr_, Arpack_jll.libarpack), Cvoid,
    (Ptr{UInt8},
     Ref{LinearAlgebra.BlasInt}, # bool
     Ref{LinearAlgebra.BlasInt}, # size
     Ptr{Float64},
     Ref{Float64}),
    whichstr, apply, n, x1, x2)
end

function arpack_dsgets(
  ishift::Int, # Input
  which::Symbol, # Input
  kev::Int, # Input
  np::Int, # Input
  ritz::Vector{Float64}, # Input/Output
  bounds::Vector{Float64}, # Input/Output
  shifts::Vector{Float64} # Input/Output
  )
  ccall((:dsgets_, Arpack_jll.libarpack), Cvoid,
    (Ref{LinearAlgebra.BlasInt},
     Ptr{UInt8},
     Ref{LinearAlgebra.BlasInt},
     Ref{LinearAlgebra.BlasInt},
     Ptr{Float64},
     Ptr{Float64},
     Ptr{Float64}),
   ishift, string(which), kev, np, ritz, bounds, shifts)
end

function arpack_dstqrb!(n::Int, d::Vector{Float64}, e::Vector{Float64},
                      z::Vector{Float64}, work::Vector{Float64})
  info = Ref{LinearAlgebra.BlasInt}(-1)
  ccall((:dstqrb_, Arpack_jll.libarpack), Cvoid,
    (Ref{LinearAlgebra.BlasInt},
     Ptr{Float64},
     Ptr{Float64},
     Ptr{Float64},
     Ptr{Float64},
     Ref{LinearAlgebra.BlasInt}),
    n, d, e, z, work, info)
  return info[]
end

function arpack_dseigt!(rnorm::T,
  n::Int, H::AbstractMatrix{T}, ldh::Int,
  eig::AbstractVecOrMat{T}, bounds::AbstractVecOrMat{T}, 
  workl::AbstractVecOrMat{T}
) where T <: Float64 
  info = Ref{LinearAlgebra.BlasInt}(-1)
  ccall((:dseigt_, Arpack_jll.libarpack), Cvoid,
    (Ref{Float64},
      Ref{LinearAlgebra.BlasInt},
      Ptr{Float64},
      Ref{LinearAlgebra.BlasInt},
      Ptr{Float64},
      Ptr{Float64},
      Ptr{Float64},
      Ref{LinearAlgebra.BlasInt}),
    rnorm, n, H, ldh, eig, bounds, workl, info
  )
  return info[]
end

##
import Arpack_jll, LinearAlgebra
function arpack_dgetv0!(ido::Ref{LinearAlgebra.BlasInt}, bmat::Symbol, itry::Int, initv::Bool,
  n::Int, j::Int, v::StridedMatrix{Float64}, ldv::Int,
    resid::StridedVecOrMat{Float64}, rnorm::Ref{Float64},
    ipntr::StridedVecOrMat{LinearAlgebra.BlasInt}, workd::StridedVecOrMat{Float64})
  ierr = Ref{LinearAlgebra.BlasInt}(0)
  # NOTE, arpack doesn't touch ierr unless ido[] == 0 or there is
  # a restart failure.
  # @show ido[], itry, initv, n, j, string(bmat)
  ccall((:dgetv0_, Arpack_jll.libarpack), Cvoid,
    (Ref{LinearAlgebra.BlasInt},
     Ptr{UInt8},
     Ref{LinearAlgebra.BlasInt},
     Ref{LinearAlgebra.BlasInt}, # really a logical...
     Ref{LinearAlgebra.BlasInt},
     Ref{LinearAlgebra.BlasInt},
     Ptr{Float64},
     Ref{LinearAlgebra.BlasInt},
     Ptr{Float64},
     Ref{Float64},
     Ptr{LinearAlgebra.BlasInt},
     Ptr{Float64},
     Ref{LinearAlgebra.BlasInt},
     Int),
    ido, string(bmat), itry, initv, n, j, v, ldv, resid, rnorm, ipntr, workd, ierr, 1)
  return ierr[]
end


##
import Arpack_jll, LinearAlgebra
function arpack_dsaitr!(
    ido::Ref{LinearAlgebra.BlasInt}, 
    bmat::Symbol, 
    n::Int, 
    k::Int,
    np::Int, 
    mode::Int, 
    resid::StridedVecOrMat{Float64},
    rnorm::Ref{Float64},
    v::StridedVecOrMat{Float64},
    ldv::Int,
    h::StridedVecOrMat{Float64},
    ldh::Int,
    ipntr::StridedVecOrMat{LinearAlgebra.BlasInt},
    workd::StridedVecOrMat{Float64}
)
  info = Ref{LinearAlgebra.BlasInt}(0)
  # NOTE, arpack doesn't touch ierr unless ido[] == 0 or there is
  # a restart failure.
  # @show ido[], itry, initv, n, j, string(bmat)
  ccall((:dsaitr_, Arpack_jll.libarpack), Cvoid,
    (Ref{LinearAlgebra.BlasInt}, # ido
     Ptr{UInt8}, # bmat
     Ref{LinearAlgebra.BlasInt}, #n
     Ref{LinearAlgebra.BlasInt}, #k 
     Ref{LinearAlgebra.BlasInt}, #np
     Ref{LinearAlgebra.BlasInt}, #mode 
     Ptr{Float64}, #resid
     Ref{Float64}, #rnorm
     Ptr{Float64}, #v
     Ref{LinearAlgebra.BlasInt}, #ldv
     Ptr{Float64}, #h
     Ref{LinearAlgebra.BlasInt}, #lhd
     Ptr{LinearAlgebra.BlasInt}, #ipntr
     Ptr{Float64}, #workd
     Ref{LinearAlgebra.BlasInt}, Int), #info
    ido, string(bmat), n, k, np, mode, resid, rnorm, 
    v, ldv, h, ldh, ipntr, workd, info, 1)
  return info[]
end


##
import Arpack_jll, LinearAlgebra
function arpack_dsaup2!(
  ido::Ref{Int}, 
  bmat::Symbol,
  n::Int,
  which::Symbol,
  nev::Ref{Int},
  np::Ref{Int}, 
  tol::Float64,
  resid::StridedVecOrMat{Float64},
  mode::Int, 
  iupd::Int,
  ishift::Int,
  mxiter::Ref{Int},
  V::StridedMatrix{Float64},
  ldv::Int, 
  H::StridedMatrix{Float64}, 
  ldh::Int,
  ritz::StridedVecOrMat{Float64},
  bounds::StridedVecOrMat{Float64},
  Q::StridedMatrix{Float64},
  ldq::Int, 
  workl::StridedVecOrMat{Float64},
  ipntr::StridedVecOrMat{Int},
  workd::StridedVecOrMat{Float64},
  info_initv0::Int, # info in Arpack, but we return info... 
)
  info = Ref{LinearAlgebra.BlasInt}(info_initv0)
  ccall((:dsaup2_, Arpack_jll.libarpack), Cvoid,
    (Ref{LinearAlgebra.BlasInt}, # ido
     Ptr{UInt8}, # bmat
     Ref{LinearAlgebra.BlasInt}, #n
     Ptr{UInt8}, # which
     Ref{LinearAlgebra.BlasInt}, # nev
     Ref{LinearAlgebra.BlasInt}, # np
     Ref{Float64}, # tol
     Ptr{Float64}, # resid
     Ref{LinearAlgebra.BlasInt}, # mode
     Ref{LinearAlgebra.BlasInt}, # iupd
     Ref{LinearAlgebra.BlasInt}, # ishift
     Ref{LinearAlgebra.BlasInt}, # mxiter
     Ptr{Float64}, # V
     Ref{LinearAlgebra.BlasInt}, # ldv
     Ptr{Float64}, # H
     Ref{LinearAlgebra.BlasInt}, # ldh
     Ptr{Float64}, # ritz
     Ptr{Float64}, # bounds 
     Ptr{Float64}, # Q
     Ref{LinearAlgebra.BlasInt}, # ldq
     Ptr{Float64}, # workl
     Ptr{LinearAlgebra.BlasInt}, # ipntr
     Ptr{Float64}, # workd
     Ref{LinearAlgebra.BlasInt},
     Int, Int), #info
    ido, string(bmat), n, string(which), nev, np, tol, 
    resid, 
    mode, iupd, ishift, mxiter, 
    V, ldv, H, ldh, ritz, bounds, Q, ldq, 
    workl, ipntr, workd, info, 1, 2)
  return info[]
end


##
import Arpack_jll, LinearAlgebra
function arpack_dsaupd!(
  ido::Ref{Int}, 
  bmat::Symbol,
  n::Int,
  which::Symbol,
  nev::Int,
  tol::Float64,
  resid::StridedVecOrMat{Float64},
  ncv::Int, 
  V::StridedMatrix{Float64},
  ldv::Int,
  iparam::StridedVecOrMat{Int},
  ipntr::StridedVecOrMat{Int},
  workd::StridedVecOrMat{Float64},
  workl::StridedVecOrMat{Float64},
  lworkl::Int,   
  info_initv0::Int, # info in Arpack, but we return info... 
)
  if tol == 0 && copysign(1.0, tol) > 0 
    @warn("This is probably a mistake because we aren't passing a ref for tol to get the arpack auto-tol\n"*
          "if this is really what you want, pass -0.0 to stop this warning")
  end 
  info = Ref{LinearAlgebra.BlasInt}(info_initv0)
  ccall((:dsaupd_, Arpack_jll.libarpack), Cvoid,
    (Ref{LinearAlgebra.BlasInt}, # ido
     Ptr{UInt8}, # bmat
     Ref{LinearAlgebra.BlasInt}, #n
     Ptr{UInt8}, # which
     Ref{LinearAlgebra.BlasInt}, # nev
     Ref{Float64}, # tol
     Ptr{Float64}, # resid
     Ref{LinearAlgebra.BlasInt}, # ncv
     Ptr{Float64}, # V
     Ref{LinearAlgebra.BlasInt}, # ldv
     Ptr{LinearAlgebra.BlasInt}, # iparam
     Ptr{LinearAlgebra.BlasInt}, # ipntr
     Ptr{Float64}, # workd
     Ptr{Float64}, # workl
     Ref{LinearAlgebra.BlasInt}, # lworkl
     Ref{LinearAlgebra.BlasInt}, # info
     Int, Int), #info
    ido, string(bmat), n, string(which), nev, tol,
    resid, ncv, 
    V, ldv, iparam, ipntr, workd, workl, lworkl, info, 1, 2)
  return info[]
end

##
import Arpack_jll, LinearAlgebra
function arpack_dsapps!(
  n::Int,
  kev::Int,
  np::Int,
  shift::StridedVecOrMat{Float64},
  V::StridedMatrix{Float64},
  ldv::Int,
  H::StridedVecOrMat{Float64},
  ldh::Int,
  resid::StridedVecOrMat{Float64},
  Q::StridedMatrix{Float64},
  ldq::Int,
  workd::StridedVecOrMat{Float64}
)
  ccall((:dsapps_, Arpack_jll.libarpack), Cvoid,
    (Ref{LinearAlgebra.BlasInt}, # n
     Ref{LinearAlgebra.BlasInt}, # kev
     Ref{LinearAlgebra.BlasInt}, # np
     Ptr{Float64}, # shift
     Ptr{Float64}, # v
     Ref{LinearAlgebra.BlasInt}, # ldv 
     Ptr{Float64}, # h
     Ref{LinearAlgebra.BlasInt}, # ldh,
     Ptr{Float64}, # resid
     Ptr{Float64}, # Q
     Ref{LinearAlgebra.BlasInt}, # ldq
     Ptr{Float64}, # workd 
     ),
    n, kev, np, shift, V, ldv, H, ldh, resid, Q, ldq, workd)
  return nothing 
end

##
import Arpack_jll, LinearAlgebra
function arpack_dseupd!(
  rvec::Bool, 
  select::Vector{Int},
  d::Vector{Float64},
  Z::Matrix{Float64},
  ldz::Int, 
  sigma::Float64, 
  bmat::Symbol,
  n::Int,
  which::Symbol,
  nev::Int,
  tol::Float64,
  resid::StridedVecOrMat{Float64},
  ncv::Int, 
  V::StridedMatrix{Float64},
  ldv::Int,
  iparam::StridedVecOrMat{Int},
  ipntr::StridedVecOrMat{Int},
  workd::StridedVecOrMat{Float64},
  workl::StridedVecOrMat{Float64},
  lworkl::Int,   
)
  info = Ref{LinearAlgebra.BlasInt}(0)
  ccall((:dseupd_, Arpack_jll.libarpack), Cvoid,
    (Ref{LinearAlgebra.BlasInt}, # rvec
     Ptr{UInt8}, # howmny, 
     Ptr{LinearAlgebra.BlasInt}, # select, 
     Ptr{Float64}, # d
     Ptr{Float64}, # Z
     Ref{LinearAlgebra.BlasInt}, # ldz
     Ref{Float64}, # sigma 
     Ptr{UInt8}, # bmat
     Ref{LinearAlgebra.BlasInt}, #n
     Ptr{UInt8}, # which
     Ref{LinearAlgebra.BlasInt}, # nev
     Ref{Float64}, # tol
     Ptr{Float64}, # resid
     Ref{LinearAlgebra.BlasInt}, # ncv
     Ptr{Float64}, # V
     Ref{LinearAlgebra.BlasInt}, # ldv
     Ptr{LinearAlgebra.BlasInt}, # iparam
     Ptr{LinearAlgebra.BlasInt}, # ipntr
     Ptr{Float64}, # workd
     Ptr{Float64}, # workl
     Ref{LinearAlgebra.BlasInt}, # lworkl
     Ref{LinearAlgebra.BlasInt}, # info
     Int, Int, Int), #info
    rvec, string(:A), select, d, Z, ldz, sigma, 
    string(bmat), n, string(which), nev, tol,
    resid, ncv, 
    V, ldv, iparam, ipntr, workd, workl, lworkl, info, 1, 1, 2)
  return info[]
end
