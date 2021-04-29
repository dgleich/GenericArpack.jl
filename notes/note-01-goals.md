The goal of this procedure is to port the double-precision ARPACK
for symmetric matrices in Julia. Including all ARPACK stuff. So this should
give "exactly" what ARPACK does but be a pure Julia implementation.
(Where exactly is ... it should be executing roughly the same sequence of
floating point operations and can differ on levels that would be expected
for different compilers compiling the same code.)

- not a goal to "Julia-ize" the package; I want to keep as close to the FORTRAN
as possible so that I might be able to replace calls to Julia's Arpack.saupd /
Arpack.seupd (which call the Fortran library) with this code.
- small internal function changes are okay, e.g. ARPACK has various debugging
and timing stuff that would need to be done differently in Julia.
- small simplifications, e.g. if a function computes a single Int, we can
rewrite that to return the Int rather than writing it into an array like in
FORTRAN.
- Why? Why not implement my own ImplicitRestart/Eigensolver/Etc.? Simple: I trust
ARPACK. Also, I want to understand exactly what the symmetric ARPACK solver is doing.

For context, I spent about a semester studying these solvers with Gene Golub
to see if there were obvious ideas to make computing the Fiedler vectors faster,
but beyond the high level algorithmics and lots of experiments, nothing came of
it. This project is different in that the end goal is easy -- port the code,
with a focus on explaining in detail what each step does.

Since it's good to have a driving use-case,
the idea is to make the following MatrixNetworks.jl code work without the
ARPACK dependency.

        """
        `_symeigs_smallest_arpack`
        ---
        Compute eigenvalues and vectors using direct calls to the ARPACK wrappers
        to get type-stability. This function works for symmetric matrices.
        This function works on Float32 and Float64-valued sparse matrices.
        It returns the smallest set of eigenvalues and vectors.
        It only works on matrices with more than 21 rows and columns.
        (Use a dense eigensolver for smaller problems.)
        Functions
        ---------
        - `(evals,evecs) = _symeigs_smallest_arpack(A::SparseMatrixCSC{V,Int},
                                nev::Int,tol::V,maxiter::Int, v0::Vector{V})`
        Inputs
        ------
        - `A`: the sparse matrix, must be symmetric
        - `nev`: the number of eigenvectors requested
        - `tol`: the relative tolerance of the eigenvalue computation
        - `maxiter`: the maximum number of restarts
        - `v0`: the initial vector for the Lanczos process
        Example
        -------
        This is an internal function.
        """
        function _symeigs_smallest_arpack(
                    A::SparseMatrixCSC{V,Int},nev::Int,tol::V,maxiter::Int,
                    v0::Vector{V}) where V

            n::Int = checksquare(A) # get the size
            @assert n >= 21

            # setup options
            mode = 1
            sym = true
            iscmplx = false
            bmat = String("I")
            ncv = min(max(2*nev,20),n-1)

            whichstr = String("SA")
            ritzvec = true
            sigma = 0.

            #TOL = Array{V}(undef,1)
            #TOL[1] = tol
            TOL = Ref(tol)
            lworkl = ncv*(ncv + 8)
            v = Array{V}(undef, n, ncv)
            workd = Array{V}(undef, 3*n)
            workl = Array{V}(undef, lworkl)
            resid = Array{V}(undef, n)

            resid[:] = v0[:]

            info = zeros(BlasInt, 1)
            info[1] = 1

            iparam = zeros(BlasInt, 11)
            ipntr = zeros(BlasInt, 11)
            ido = zeros(BlasInt, 1)

            iparam[1] = BlasInt(1)
            iparam[3] = BlasInt(maxiter)
            iparam[7] = BlasInt(mode)

            # this is a helpful indexing vector
            zernm1 = 0:(n-1)

            while true
                # This is the reverse communication step that ARPACK does
                # we need to extract the desired vector and multiply it by A
                # unless the code says to stop
                Arpack.saupd(
                    ido, bmat, n, whichstr, nev, TOL, resid, ncv, v, n,
                    iparam, ipntr, workd, workl, lworkl, info)

                load_idx = ipntr[1] .+ zernm1
                store_idx = ipntr[2] .+ zernm1

                x = workd[load_idx]

                if ido[1] == 1
                    workd[store_idx] = A*x
                elseif ido[1] == 99
                    break
                else
                    error("unexpected ARPACK behavior")
                end
            end

            # Once the process terminates, we need to extract the
            # eigenvectors.

            # calls to eupd
            howmny = String("A")
            select = Array{BlasInt}(undef, ncv)

            d = Array{V}(undef, nev)
            sigmar = ones(V,1)*sigma
            ldv = n
            Arpack.seupd(ritzvec, howmny, select, d, v, ldv, sigmar,
                bmat, n, whichstr, nev, TOL, resid, ncv, v, ldv,
                iparam, ipntr, workd, workl, lworkl, info)
            if info[1] != 0
                error("unexpected ARPACK exception")
            end

            # Now we want to return them in sorted order (smallest first)
            p = sortperm(d)

            d = d[p]
            vectors = v[1:n,p]

            return (d,vectors,iparam[5])
        end
