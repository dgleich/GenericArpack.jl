Me on learning about `save` variables in Fortran and accessing them from Julia. (Part ... something ... in my ongoing segment on porting ARPACK to Julia)
=========

Huh, what are these things in a FORTRAN code marked `save` ...
[Oh... jeez, is that persistent state in a subroutine?](https://stackoverflow.com/questions/2893097/fortran-save-statement) i.e. like global variables? ugh, that is going to be annoying to port. (More on a strategy in a bit.)

e.g. [dgetv0.f in ARPACK](https://github.com/opencollab/arpack-ng/blob/master/SRC/dgetv0.f)

      save       first, iseed, inits, iter, msglvl, orth, rnorm0

How could these possibly be implemented in a shared library? Well, most likely they are just fixed offsets in memory somewhere that the subroutine knows about and manipulates. So can I see them from Julia?

- Try dumping the dylib/so file to find if these are listed in there. (They ought to be as they should take up persistent locations in memory...) (Any of the following commands show them...)

        nm -a libarpack.2.0.0.dylib
        objdump -p libarpack.2.0.0.dylib
        dsymutil -s libarpack.2.0.0.dylib
        nm -nm libarpack.2.0.0.dylib

- They are! e.g. In this case, I picked `iseed` which is the ARPACK random number seed/state.

        dgleich@circulant test-arpack-state % nm -nm /Users/dgleich/.julia/artifacts/e1631d3aec690feb8f2330c629e853190df49fbe/lib/libarpack.2.1.0.dylib | grep iseed
        000000000005d060 (__DATA,__bss5) non-external _iseed.3637
        000000000005d080 (__DATA,__bss5) non-external _iseed.3636
        000000000005d0a0 (__DATA,__bss5) non-external _iseed.3636
        000000000005d0c0 (__DATA,__bss5) non-external _iseed.3637

   Well, there are a few things listed here... but probably something
   we can resolve via code.

- Try and lookup symbols via `dlsym`...

        using Arpack
        using Arpack_jll
        using LinearAlgebra
        libar = Base.Libc.dlopen(Arpack_jll.libarpack)
        ## Try and get offset of private symbol directly
        offset = Base.Libc.dlsym(libar, "iseed.3637")
        offset = Base.Libc.dlsym(libar, "_iseed.3637")
        offset = Base.Libc.dlsym(libar, "iseed")

    doesn't work because `dlsym only
    has exported/public symbols accessible. (Why, I don't know because they are in the darn file, so it's like they want to make it harder to access things... which I guess is a good thing?)

- Try and look at offsets from the load base. e.g. if the dylib/so file says that function `dgetv0` is at address 0x0010f00 and I load that at address 0x5020e00, then the offset to other addresses should be the same... which it is for functions (easy to verify)

        dgleich@circulant test-arpack-state % nm -nm /Users/dgleich/.julia/artifacts/e1631d3aec690feb8f2330c629e853190df49fbe/lib/libarpack.2.1.0.dylib | grep getv0
        nm -nm /Users/dgleich/.julia/artifacts/e1631d3aec690feb8f2330c629e853190df49fbe/lib/libarpack.2.1.0.dylib | grep getv0

        000000000001f2a0 (__TEXT,__text) external _dgetv0_

        julia> offset = Base.Libc.dlsym(libar, "dgetv0_")
        Ptr{Nothing} @0x0000000170bac2a0

        julia> base_offset = offset - 0x000000000001f2a0
        Ptr{Nothing} @0x0000000170b8d000

   So we have a base offset, which means that other functions and symbols should respect that...

        dgleich@circulant test-arpack-state % nm -nm /Users/dgleich/.julia/artifacts/e1631d3aec690feb8f2330c629e853190df49fbe/lib/libarpack.2.1.0.dylib | grep arscnd
        0000000000008a00 (__TEXT,__text) external _arscnd_

        julia> arscnd_offset = Base.Libc.dlsym(libar, "arscnd_")
        Ptr{Nothing} @0x0000000170b95a00

        julia> arscnd_offset == 0x0000000000008a00 + base_offset
        true

    Awesome, so let's check out some of those iseed variables...

        julia> iseed_offsets = [0x5d060, 0x5d080, 0x5d0a0, 0x5d0c0]
        julia> iseeds = Ptr{Int64}.(base_offset .+ iseed_offsets)

        julia> unsafe_load.(iseeds)
        4-element Vector{Int64}:
         0
         0
         0
         0

        julia> unsafe_load.(iseeds .+ 8)
        4-element Vector{Int64}:
         0
         0
         0
         0

    Okay, so these iseed variables are currently 0, because we haven't
    run anything yet. (Internally in `dgetv0` they are used to track
    the state of the random number generator... )

        julia> using Arpack; eigs(randn(50,50))
        julia> unsafe_load.(iseeds)

    Hmm. That shows the same values. Maybe another set? ARPACK also uses
    `first` as some one-time routine initialization tracked with `save`
    variables.

        $ nm -nm /Users/dgleich/.julia/artifacts/e1631d3aec690feb8f2330c629e853190df49fbe/lib/libarpack.2.1.0.dylib | grep first

        julia> first_offsets="""000000000005c398 (__DATA,__data) non-external _first.3636
        000000000005c3a0 (__DATA,__data) non-external _first.3637
        000000000005c3b0 (__DATA,__data) non-external _first.3635
        000000000005c3b8 (__DATA,__data) non-external _first.3640
        000000000005c3c0 (__DATA,__data) non-external _first.3634
        000000000005c3c8 (__DATA,__data) non-external _first.3639
        000000000005c3d8 (__DATA,__data) non-external _first.3635
        000000000005c3e0 (__DATA,__data) non-external _first.3640
        000000000005c3e8 (__DATA,__data) non-external _first.3634
        000000000005c3f0 (__DATA,__data) non-external _first.3639
        000000000005c400 (__DATA,__data) non-external _first.3636
        000000000005c408 (__DATA,__data) non-external _first.3637
        000000000005c780 (__DATA,__bss3) non-external _first.3633
        000000000005c8d8 (__DATA,__bss3) non-external _first.3632
        000000000005cc08 (__DATA,__bss3) non-external _first.3632
        000000000005cea8 (__DATA,__bss3) non-external _first.3633""" |>
          x->split(x,"\n") |>
          x->map(y->y[1:16],x) |>
          x->parse.(UInt64, x[1:16];base=16)

    (Sorry for syntax fancyness, `"0x001" |> x->parse(Int,x)` just creates a little
    mini-function and chains a set of function evaluations. Handy for little
    sequences like that where I split, then break apart, then parse each.)

        julia> firsts = Ptr{Int64}.(base_offset .+ first_offsets)
        julia> unsafe_load.(firsts)
        16-element Vector{Int64}:
         1
         1
         1
         1
         1
         1
         1
         1
         1
         1
         1
         1
         0
         0
         0
         0

    Which means that everything hasn't really been run yet. (This is the
    same right after loading the library. The 0's at the end
    are there as well.)

- So, hmm... that doesn't seem to work. (Or at least I don't see
  changes I expect to see...) But still the function addresses are exactly
  where I predict them.

- Try and check out dissambled routines for help. [Learn how to get dissassembled
  code from `objdump`](https://stackoverflow.com/questions/22769246/how-to-disassemble-one-single-function-using-objdump)

        dgleich@circulant test-arpack-state % objdump -D /Users/dgleich/.julia/artifacts/e1631d3aec690feb8f2330c629e853190df49fbe/lib/libarpack.2.1.0.dylib | less

        000000000001f2a0 _dgetv0_:
           1f2a0: 41 57                         pushq   %r15
           1f2a2: 41 56                         pushq   %r14
           1f2a4: 49 89 d7                      movq    %rdx, %r15
           1f2a7: 41 55                         pushq   %r13
           1f2a9: 41 54                         pushq   %r12
           1f2ab: 49 89 fd                      movq    %rdi, %r13
           1f2ae: 55                            pushq   %rbp
           1f2af: 53                            pushq   %rbx
           1f2b0: 49 89 f4                      movq    %rsi, %r12
           1f2b3: 4c 89 c3                      movq    %r8, %rbx
           1f2b6: 4d 89 ce                      movq    %r9, %r14
           1f2b9: 48 83 ec 28                   subq    $40, %rsp
           1f2bd: 48 83 3d e3 d0 03 00 00       cmpq    $0, 250083(%rip)
           1f2c5: 74 29                         je      41 <_dgetv0_+0x50>
           1f2c7: 66 0f 6f 05 31 f4 02 00       movdqa  193585(%rip), %xmm0
           1f2cf: 48 c7 05 ce d0 03 00 00 00 00 00      movq    $0, 250062(%rip)
           1f2da: 0f 29 05 9f dd 03 00          movaps  %xmm0, 253343(%rip)
           1f2e1: 66 0f 6f 05 27 f4 02 00       movdqa  193575(%rip), %xmm0
           1f2e9: 0f 29 05 a0 dd 03 00          movaps  %xmm0, 253344(%rip)
           1f2f0: 49 83 7d 00 00                cmpq    $0, (%r13)
           1f2f5: 0f 84 85 02 00 00             je      645 <_dgetv0_+0x2e0>
           1f2fb: 41 0f b6 04 24                movzbl  (%r12), %eax

    In case you want a spoiler, it's  `cmpq    $0, 250083(%rip)` that
    checks "if first == 1" on the first subroutine call
    and so the address of first is `0x1f2c5 + 250083`
    (the NEXT instruction offset + `%rip offset`) so the offset is actually
    `0x0005c3a8` (just add the numbers and write in hex).

    (This is actually off by just a bit because of how Fortran seems
    to store logicals. There are a few details I don't fully understand
    but it's close enough I'm not worried. )

- Next up. Check code on Linux, the julia code I worked on above
  seems to work there????!??!? Everything I was doing above
  should really have worked?? Amazing. (Just had to pull new offsets
  from that `libarpack`... )

- But wait, gosh,  why doesn't this work on the mac I was sitting in front of?
  Is this some weirdo Mac dlsym loader/protection system that
  randomizes offsets... and keeps it from working on the mac?

- Develop C/C++ code to do the same, verify that the addresses I'm getting
  for functions seem correct. (Meanwhile, learn how to use `lldb`...!)
  Okay, everything `lldb` seems to be doing is exactly what I think ought to
  work. (See the C/C++ code below.)

- Check that I can use the C++ code and pointer offsets as I'm doing in
  Julia and on the mac in order to

- [Learn more details about position independent code and why everything
is in terms of offsets from the instruction pointer `%rip`.](https://eli.thegreenplace.net/2011/11/11/position-independent-code-pic-in-shared-libraries-on-x64) (See that disassembled code above
    for an example)

- Start learning about ccall, julia library loading. (Hoping to find some pointer
    as to why there would be this weirdness... )

- [Run into this thread on `dlsym`,`dlopen`](https://github.com/JuliaLang/julia/issues/23459)

- Double check to see that I've gotten the right library loaded...

        julia> Base.Libc.dllist()
        152-element Vector{String}:
         "/usr/lib/libSystem.B.dylib"
         "/Applications/Julia-1.6.app/Contents/Resources/julia/lib/libjulia.1.6.dylib"
         "/usr/lib/system/libcache.dylib"
         "/usr/lib/system/libcommonCrypto.dylib"
         "/usr/lib/system/libcompiler_rt.dylib"
         "/usr/lib/system/libcopyfile.dylib"
         "/usr/lib/system/libcorecrypto.dylib"
         "/usr/lib/system/libdispatch.dylib"
         "/usr/lib/system/libdyld.dylib"
         "/usr/lib/system/libkeymgr.dylib"
         ⋮
         "/Applications/Julia-1.6.app/Contents/Resources/julia/lib/julia/libmbedx509.2.24.0.dylib"
         "/Applications/Julia-1.6.app/Contents/Resources/julia/lib/julia/libmbedcrypto.2.24.0.dylib"
         "/Applications/Julia-1.6.app/Contents/Resources/julia/lib/julia/libdSFMT.dylib"
         "/System/Library/PrivateFrameworks/DebugSymbols.framework/Versions/A/DebugSymbols"
         "/System/Library/PrivateFrameworks/CoreServicesInternal.framework/Versions/A/CoreServicesInternal"
         "/Applications/Julia-1.6.app/Contents/Resources/julia/lib/julia/libstdc++.6.dylib"
         "/Applications/Julia-1.6.app/Contents/Resources/julia/lib/julia/libgomp.1.dylib"
         "/Users/dgleich/.julia/artifacts/e1631d3aec690feb8f2330c629e853190df49fbe/lib/libarpack.2.1.0.dylib"
         "/Users/dgleich/.julia/packages/Arpack/zCmTA/deps/usr/lib/libarpack.2.0.0.dylib"

- Oh...

- You see that? In the bottom there? There are _two_ `libarpack`s listed.
- I was looking at the wrong one. :(
- 🤦

- Just for kicks now...

        dgleich@circulant test-arpack-state % nm -nm /Users/dgleich/.julia/packages/Arpack/zCmTA/deps/usr/lib/libarpack.2.0.0.dylib | grep iseed
        000000000005afa0 (__DATA,__bss5) non-external _iseed.3637
        000000000005afc0 (__DATA,__bss5) non-external _iseed.3636
        000000000005afe0 (__DATA,__bss5) non-external _iseed.3636
        000000000005b000 (__DATA,__bss5) non-external _iseed.3637

        julia> libar = Base.Libc.dlopen("/Users/dgleich/.julia/packages/Arpack/zCmTA/deps/usr/lib/libarpack.2.0.0.dylib")
        julia> offset = Base.Libc.dlsym(libar, "dgetv0_")
        Ptr{Nothing} @0x000000017596f730

        dgleich@circulant test-arpack-state % nm -nm /Users/dgleich/.julia/packages/Arpack/zCmTA/deps/usr/lib/libarpack.2.0.0.dylib | grep dgetv0
        000000000001e730 (__TEXT,__text) external _dgetv0_

        julia> base_offset = offset - 0x000000000001e730
        Ptr{Nothing} @0x0000000175951000

        julia> libar = Base.Libc.dlopen(Arpack_jll.libarpack)
        julia> iseed_offsets = [0x5afa0, 0x5afc0, 0x5afe0, 0x5b000]
        julia> iseeds = Ptr{Int64}.(base_offset .+ iseed_offsets)

        julia> unsafe_load.(iseeds)

        julia> unsafe_load.(iseeds)
        4-element Vector{Int64}:
            0
         1626
            0
            0

        julia> eigs(randn(50,50))
        julia>  unsafe_load.(iseeds)
        4-element Vector{Int64}:
           0
         561
           0
           0

- 🤦 indeed. (This took way more time than I would have liked. Approx. a full day.)

- Should have trusted the fact that this all worked in Linux, where somehow `Arpack` is using the same library as `Arpack_jll`

- One other thing I learned. On the macos build of arpack. There is no
timing infomation collected. The `arscnd` function (which should return
some time value in seconds) just is a no-op.

        0000000000008a00 _arscnd_:
            8a00: c7 07 00 00 00 00             movl    $0, (%rdi)
            8a06: c3                            retq
            8a07: 90                            nop
            8a08: 90                            nop
            8a09: 90                            nop
            8a0a: 90                            nop
            8a0b: 90                            nop
            8a0c: 90                            nop
            8a0d: 90                            nop
            8a0e: 90                            nop
            8a0f: 90                            nop

    (Or rather, this seems to write 0 to the location passed
    into `arscnd`!)

    I learned this, of course, because one of the thing `arpack` does
    is collect timing information to persist between function calls
    so it can time user-code via the call-backs and such.

    So in debugging, I was often looking at locations where time
    information ought to be stored. Finding it continually
    empty was ... surprising.

- Another thing I learned, you can access the common blocks of FORTRAN code
  via `cglobal` in Julia.

        julia> debug = cglobal(("debug_", "/Users/dgleich/.julia/packages/Arpack/zCmTA/deps/usr/lib/libarpack.2.0.0.dylib"),Int64)
        Ptr{Int64} @0x00000001759ac020

    [Documentation and structure of debug block here.](https://github.com/opencollab/arpack-ng/blob/master/SRC/debug.h)
    But we can use it to turn on lots of fun ARPACK debug information!

        julia> unsafe_store!(debug, 2, 11) # mnaupd turn on debugging for mnaupd (non-symmetric solver)
        julia> unsafe_store!(debug, -14, 2) # ndigit - use 14 digits of precision (ndigit)
        julia> unsafe_store!(debug, 4, 3) # mgetv0 show lots of information from dgetv0
        julia> unsafe_store!(debug, 6, 1) # logfile set logfil to 6, the default stdout

        julia> eigs(randn(50,50))

         _getv0: B-norm of initial / restarted starting vector
         -----------------------------------------------------
            1 -    1:      4.1596362767673D+00


         _getv0: initial / restarted starting vector
         -------------------------------------------
            1 -    2:     -8.3153140665625D-01     8.3995344740674D-01
            3 -    4:      9.6704039685513D-01     3.6866117348409D-01
            5 -    6:      6.0284542702433D-03     2.0318894267245D-01
            7 -    8:     -1.3495066933914D-01     4.1078102410386D-01
            9 -   10:      4.3888096465232D-01    -5.7522033107717D-01
           11 -   12:     -8.7293721104508D-01     9.8788590316195D-01
           13 -   14:     -1.2282299641820D-01     7.4369927039981D-01
           15 -   16:      8.7906957254355D-01    -4.7762629697247D-01
           17 -   18:     -4.7162805784843D-02     1.1706962983657D-01
           19 -   20:      7.8492285203887D-01     2.6127819008416D-01
           21 -   22:     -8.2777335418346D-01    -6.5745389045117D-01
           23 -   24:      4.4467420164616D-01     9.9532171436879D-01
           25 -   26:      6.9740537844088D-01    -9.2167903292103D-01
           27 -   28:     -9.1859980666224D-01    -7.1191753019861D-01
           29 -   30:     -1.5899117215995D-01     1.5815438448009D-01
           31 -   32:     -5.1806821604240D-01     3.8970905952163D-01
           33 -   34:     -8.7740030659049D-01    -5.7437967251424D-02
           35 -   36:      7.9785742115958D-01     3.1842249389620D-01
           37 -   38:     -6.1567371152282D-01    -2.1826976841397D-01
           39 -   40:      8.0618226869270D-01     2.7929452805704D-01
           41 -   42:     -7.8813902002933D-01     1.7058205959301D-01
           43 -   44:     -4.0247429259527D-01    -6.0433106862224D-01
           45 -   46:     -4.1118723651662D-01    -3.4363227883103D-01
           47 -   48:     -4.4783519455385D-01     2.1286373092316D-02
           49 -   50:      2.9004800999430D-01     4.4388758381478D-01


         _naupd: Number of update iterations taken
         -----------------------------------------
            1 -    1:              17


         _naupd: Number of wanted "converged" Ritz values
         ------------------------------------------------
            1 -    1:               6


         _naupd: Real part of the final Ritz values
         ------------------------------------------
            1 -    2:     -5.7794067256680D+00    -5.7794067256680D+00
            3 -    4:     -7.2262688525223D+00     6.8333693387731D+00
            5 -    6:      6.8333693387731D+00     7.6197283464233D+00


         _naupd: Imaginary part of the final Ritz values
         -----------------------------------------------
            1 -    2:      2.7424220041744D+00    -2.7424220041744D+00
            3 -    4:      0.0000000000000D+00     2.3888574879058D+00
            5 -    6:     -2.3888574879058D+00     0.0000000000000D+00


         _naupd: Associated Ritz estimates
         ---------------------------------
            1 -    2:      5.0738750651244D-16     5.0738750651244D-16
            3 -    4:      2.5131234631038D-25     1.2315725733964D-30
            5 -    6:      1.2315725733964D-30     1.0152776765181D-35



             =============================================
             = Nonsymmetric implicit Arnoldi update code =
             = Version Number:  2.4                      =
             = Version Date:    07/31/96                 =
             =============================================
             = Summary of timing statistics              =
             =============================================


             Total number update iterations             =    17
             Total number of OP*x operations            =   206
             Total number of B*x operations             =     0
             Total number of reorthogonalization steps  =    35
             Total number of iterative refinement steps =     0
             Total number of restart steps              =     0
             Total time in user OP*x operation          =     0.000000
             Total time in user B*x operation           =     0.000000
             Total time in Arnoldi update routine       =     0.000000
             Total time in naup2 routine                =     0.000000
             Total time in basic Arnoldi iteration loop =     0.000000
             Total time in reorthogonalization phase    =     0.000000
             Total time in (re)start vector generation  =     0.000000
             Total time in Hessenberg eig. subproblem   =     0.000000
             Total time in getting the shifts           =     0.000000
             Total time in applying the shifts          =     0.000000
             Total time in convergence testing          =     0.000000
             Total time in computing final Ritz vectors =     0.000000

Other references
- <https://github.com/conan-io/conan-center-index/issues/4696>
- <https://stackoverflow.com/questions/22769246/how-to-disassemble-one-single-function-using-objdump>
- <https://eli.thegreenplace.net/2011/11/11/position-independent-code-pic-in-shared-libraries-on-x64>
- <https://github.com/conan-io/conan-center-index/issues/4696>


C/C++ Code to check
-------------------

Compile with (edit to your own paths)

    clang -g -fPIC test-arpack-state-dlopen.c \
    -Wl,-rpath /Users/dgleich/.julia/artifacts/e1631d3aec690feb8f2330c629e853190df49fbe/lib \
    -Wl,-rpath /Applications/Julia-1.6.app/Contents/Resources/julia/lib/julia \
        -o mytest-dlopen

And here's the code. This just calls ARPACK once to see what happens.
In this case, I figured out exactly which offset to `first` to use (it
    was the first one, hah!)

    /* test_arpack_state */
    #include <inttypes.h>
    typedef int64_t blasint;
    typedef void (*dsaupd_ptr)(
      blasint* ido,
      char* bmat,
      blasint* n,
      char* whichstr,
      blasint* nev,
      double* tol,
      double* resid,
      blasint* ncv,
      double* v,
      blasint* ldv,
      blasint* iparam,
      blasint* ipntr,
      double* workd,
      double* workl,
      blasint* lworkl,
      blasint* info);

    #include <printf.h>
    #include <dlfcn.h>

    int main(int argc, char** argv) {

      blasint nev = 3;
      double tol = -1.0;
      blasint maxiter = 500;

      blasint n = 25;
      double v0[25] = {0};

      blasint ncv = 20;
      double v[25*20];
      double workd[25*3];
      double workl[20*(20+8)];
      double resid[25];
      blasint lworkl = 20*(20+8);

      blasint info = 0;

      blasint iparam[11] = {0};
      blasint ipntr[11] = {0};
      blasint ido = 0;

      iparam[0] = 1;
      iparam[2] = maxiter;
      iparam[6] = 1;

      char bmat[] = "I";
      char whichstr[] = "SA";

      blasint ldv = 25;

      void* libar = dlopen("/Users/dgleich/.julia/artifacts/e1631d3aec690feb8f2330c629e853190df49fbe/lib/libarpack.2.1.0.dylib", RTLD_LAZY|RTLD_GLOBAL);

      if (libar) {
        printf("Loaded libarpack\n");
      }

      dsaupd_ptr pfun = dlsym(libar, "dsaupd_");

      if (pfun) {
        printf("Loaded dsaupd_\n");
      }

      void* base_offset = pfun - 0x000000000002acb0; // offset from objdump.
      double* dsaupd_t0 = base_offset + 0x000000000005c604;  // offset from objdump.
      blasint* logfil = base_offset + 0x000000000005c420;
      blasint* ndig = base_offset + 0x000000000005c420+8;
      blasint* mgetv0 = base_offset + 0x000000000005c420+16;
      blasint* msaupd = base_offset + 0x000000000005c420+24;
      *mgetv0 = 2;
      *msaupd = 2;
      blasint* dgetv0_first1 = base_offset + 0x000000000005c3a0+8;
      blasint* dgetv0_first2 = base_offset + 0x000000000005c408+8;
      printf("logfil = %i\n", (int)*logfil);
      printf("ndigit = %i\n", (int)*ndig);
      printf("dgetv0_first1 = %i\n", (int)*dgetv0_first1);
      printf("dgetv0_first2 = %i\n", (int)*dgetv0_first2);

      pfun(&ido, bmat, &n, whichstr, &nev, &tol, resid, &ncv, v, &ldv, iparam,
        ipntr, workd, workl, &lworkl, &info );

      printf("base_offset = %p\n", base_offset);
      printf("dsapud = %p\n", pfun);
      printf("dsapud_t0 = %p\n", dsaupd_t0);
      printf("info = %i\n", (int)info);
      printf("ido = %i\n", (int)ido);
      printf("dsaupd_t0 = %lf\n", (double)*dsaupd_t0);

      printf("dgetv0_first1 = %i\n", (int)*dgetv0_first1);
      printf("dgetv0_first2 = %i\n", (int)*dgetv0_first2);

      //dlclose(libar);

      return 0;
    }
