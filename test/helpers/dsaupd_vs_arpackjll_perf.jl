## test allocations
using BenchmarkTools
using GenericArpack
using LinearAlgebra
using Arpack_jll

import Arpack_jll, LinearAlgebra
import SHA

LinearAlgebra.BLAS.set_num_threads(1)

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

  else
    errmsg = """Unknown libarpack sha $(bytes2hex(sha)) for $(Arpack_jll.libarpack_path)

    Please post a new issue on the github page for GenericArpack."""
    @error(errmsg)
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

function arpack_dsaupd!(
  ido::Ref{Int}, 
  bmat::String,
  n::Int,
  which::String,
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
    ido, (bmat), n, (which), nev, tol,
    resid, ncv, 
    V, ldv, iparam, ipntr, workd, workl, lworkl, info, 1, 2)
  return info[]
end

function eigrun(op,ido, ::Val{BMAT}, n, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr, workd, workl, lworkl, info_initv, state; idonow=false)  where BMAT
  niter = 0 
  nbytes = 0 

  if idonow # just try everything in one go! 
    nbytes += @allocated GenericArpack.dsaupd!(ido, Val(BMAT), n, which, nev, tol, resid, ncv, V, ldv, iparam,
      ipntr, workd, workl, lworkl, info_initv;
      state, idonow = op
    )
    if ido[] != 99
      @warn("eigrun with idonow gave unfinished result")
    end
    return niter, nbytes
  end

  while ido[] != 99
    nbytes += @allocated GenericArpack.dsaupd!(ido, Val(BMAT), n, which, nev, tol, resid, ncv, V, ldv, iparam,
      ipntr, workd, workl, lworkl, info_initv;
      state 
    )
    if ido[] == 1 || ido[] == -1
      niter += 1
      GenericArpack._i_do_now_opx_1!(op, ipntr, workd, n)
    elseif ido[] == 99
      break
    else
      @error("this only supports standard eigenvalue problems")
    end 
  end
  return niter, nbytes
end 

function eigrun_arpackjll(op, ido, bmat, n, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr, workd, workl, lworkl, info_initv)
  niter = 0 
  nbytes = 0 
  while ido[] != 99
    nbytes += @allocated arpack_dsaupd!(ido, bmat, n, which, nev, tol, resid, ncv, V, ldv, iparam,
      ipntr, workd, workl, lworkl, info_initv
    )
    
    info_initv = 0 # important for arpack... 
    if ido[] == 1 || ido[] == -1
      niter += 1
      GenericArpack._i_do_now_opx_1!(op, ipntr, workd, n)
    elseif ido[] == 99
      break
    else
      @error("this only supports standard eigenvalue problems")
    end 
  end
  return niter, nbytes
end 


##
println("Profiling arpack with reverse communcation")
begin 
  @btime begin
    eigrun(op, ido, bmat, n, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr, workd, workl, lworkl, info_initv, state);
  end setup=begin
    op = GenericArpack.ArpackSimpleOp(Diagonal(1.0:10^3))
    nev = 6
    ido = Ref{Int}(0)
    bmat = Val(:I)
    n = size(op.A,1)
    which = :LM
    tol = eps(Float64)/2 # just use the default
    resid = zeros(n)
    ncv = min(2nev, n-1)
    V = zeros(n,ncv)
    ldv = n
    mode = 1 
    iparam = zeros(Int,11)
    iparam[1] = 1
    #iparam[3] = 300 # max iteration
    iparam[3] = 300 # 
    iparam[4] = 1
    iparam[7] = mode 
    ipntr = zeros(Int,11)
    workd = zeros(3n)
    lworkl = ncv*ncv + 8*ncv
    workl = zeros(lworkl)

    info_initv = 0

    # Note that we cannot run two sequences at once and check them where we start a whole
    # second arpack call because of the expected Arpack state. 
    state = GenericArpack.ArpackState{Float64}()

    niter = 0 
  end
end

##
println("Profiling arpack with idonow ops")
begin 
  @btime begin
    eigrun(op, ido, bmat, n, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr, workd, workl, lworkl, info_initv, state; idonow=true);
  end setup=begin
    op = GenericArpack.ArpackSimpleOp(Diagonal(1.0:10^3))
    nev = 6
    ido = Ref{Int}(0)
    bmat = Val(:I)
    n = size(op.A,1)
    which = :LM
    tol = eps(Float64)/2 # just use the default
    resid = zeros(n)
    ncv = min(2nev, n-1)
    V = zeros(n,ncv)
    ldv = n
    mode = 1 
    iparam = zeros(Int,11)
    iparam[1] = 1
    #iparam[3] = 300 # max iteration
    iparam[3] = 300 # 
    iparam[4] = 1
    iparam[7] = mode 
    ipntr = zeros(Int,11)
    workd = zeros(3n)
    lworkl = ncv*ncv + 8*ncv
    workl = zeros(lworkl)

    info_initv = 0

    # Note that we cannot run two sequences at once and check them where we start a whole
    # second arpack call because of the expected Arpack state. 
    state = GenericArpack.ArpackState{Float64}()

    niter = 0 
  end
end

##
println("Profiling arpack with Arpack_jll")
begin 
  @btime begin
    eigrun_arpackjll(op, ido, bmat, n, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr, workd, workl, lworkl, info_initv);
  end setup=begin
    op = GenericArpack.ArpackSimpleOp(Diagonal(1.0:10^3))
    nev = 6
    ido = Ref{Int}(0)
    bmat = string(:I)
    n = size(op.A,1)
    which = string(:LM)
    tol = eps(Float64)/2 # just use the default
    resid = zeros(n)
    ncv = min(2nev, n-1)
    V = zeros(n,ncv)
    ldv = n
    mode = 1 
    iparam = zeros(Int,11)
    iparam[1] = 1
    iparam[3] = 300 # 
    iparam[4] = 1
    iparam[7] = mode 
    ipntr = zeros(Int,11)
    workd = zeros(3n)
    lworkl = ncv*ncv + 8*ncv
    workl = zeros(lworkl)
    info_initv = 0

    _reset_libarpack_dgetv0_iseed()
  end
end

## profiling 
using Revise, GenericArpack, LinearAlgebra, BenchmarkTools
begin
  op = GenericArpack.ArpackSimpleOp(Diagonal(1.0:10^3))
  nev = 6
  ido = Ref{Int}(0)
  bmat = Val(:I)
  n = size(op.A,1)
  which = :LM
  tol = eps(Float64)/2 # just use the default
  resid = zeros(n)
  ncv = min(2nev, n-1)
  V = zeros(n,ncv)
  ldv = n
  mode = 1 
  iparam = zeros(Int,11)
  iparam[1] = 1
  #iparam[3] = 300 # max iteration
  iparam[3] = 300 # 
  iparam[4] = 1
  iparam[7] = mode 
  ipntr = zeros(Int,11)
  workd = zeros(3n)
  lworkl = ncv*ncv + 8*ncv
  workl = zeros(lworkl)

  info_initv = 0

  # Note that we cannot run two sequences at once and check them where we start a whole
  # second arpack call because of the expected Arpack state. 
  state = GenericArpack.ArpackState{Float64}()

  niter = 0 
  @profview eigrun(op, ido, bmat, n, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr, workd, workl, lworkl, info_initv, state)
end