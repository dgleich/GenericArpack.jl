@testset "bad cases" begin
  using SparseArrays, LinearAlgebra
  function loadmat(colptr, rowval, nzval)
    n = length(colptr)-1
    SparseMatrixCSC(n,n,colptr, rowval, nzval)
  end
  A = loadmat([1, 9, 16, 22, 27, 33, 38, 47, 51, 58, 62, 71, 78, 85, 89, 95], [2, 3, 4, 7, 8, 9, 11, 15, 1, 3, 5, 7, 11, 12, 14, 1, 2, 4, 5, 6, 9, 1, 3, 4, 6, 13, 2, 3, 6, 7, 12, 15, 3, 4, 5, 10, 13, 1, 2, 5, 8, 9, 10, 11, 12, 13, 1, 7, 9, 11, 1, 3, 7, 8, 13, 14, 15, 6, 7, 11, 15, 1, 2, 7, 8, 10, 12, 13, 14, 15, 2, 5, 7, 11, 12, 13, 14, 4, 6, 7, 9, 11, 12, 15, 2, 9, 11, 12, 1, 5, 9, 10, 11, 13], [1.0619322548149914, 0.8766879002200858, -0.9064189303728223, 0.21600317366983124, -1.283210593569618, 0.94558609230275, 1.9155428985937315, 0.9974159980590636, 1.0619322548149914, 0.07774792930492884, 0.2649798493656264, 0.6491820935668625, -0.6188565210918642, 0.33167352440104847, 0.06581708949429992, 0.8766879002200858, 0.07774792930492884, -0.831030685702723, -0.9593728263464558, -1.0783982816827704, 1.5537090724942313, -0.9064189303728223, -0.831030685702723, -2.73568483736268, -0.2601662488800133, -2.0138536457007397, 0.2649798493656264, -0.9593728263464558, 1.2063734934686803, -0.23345707234182503, -0.598980265258237, -3.235531538166721, -1.0783982816827704, -0.2601662488800133, 1.2063734934686803, -1.2968161204379096, -0.5037853135483386, 0.21600317366983124, 0.6491820935668625, -0.23345707234182503, -1.635181006451923, 1.2505504541996117, -1.8979584955867908, 0.8078340897940672, 0.3488529638883331, -0.8853907377562747, -1.283210593569618, -1.635181006451923, -1.0402874079575923, -0.023802716431538706, 0.94558609230275, 1.5537090724942313, 1.2505504541996117, -1.0402874079575923, -0.4150736664751938, -0.7608617919970803, -0.1985227463733881, -1.2968161204379096, -1.8979584955867908, 0.4120271156295672, -0.033303822949929375, 1.9155428985937315, -0.6188565210918642, 0.8078340897940672, -0.023802716431538706, 0.4120271156295672, -1.7060217210513384, 0.6923442199350044, 0.5061616454994609, 0.5200080839628756, 0.33167352440104847, -0.598980265258237, 0.3488529638883331, -1.7060217210513384, -0.4259291978517505, -0.7347260397652922, -0.45812838793146876, -2.0138536457007397, -0.5037853135483386, -0.8853907377562747, -0.4150736664751938, 0.6923442199350044, -0.7347260397652922, -0.28756151077851555, 0.06581708949429992, -0.7608617919970803, 0.5061616454994609, -0.45812838793146876, 0.9974159980590636, -3.235531538166721, -0.1985227463733881, -0.033303822949929375, 0.5200080839628756, -0.28756151077851555])
  k = 3
  Avals,Avecs = eigen(Matrix(A))
  vals, vecs = eigs(Symmetric(A), k; debug=((ArpackDebug())))
  @test vals ≈ sort(sort(Avals, by=abs,rev=true)[1:k])
    
end

@testset "interface" begin 
  using SparseArrays, LinearAlgebra, StableRNGs
  n = 15 
  rng = StableRNG(1234)
  A = sprandn(rng, n, n, 0.25) |> A -> A + A' 
  Avals,Avecs = eigen(Matrix(A))
  #@show A.colptr, A.rowval, A.nzval

  k = 3 
  #vals, vecs = eigs(Symmetric(A), k; debug=set_debug_high!(ArpackDebug()))
  vals, vecs = eigs(Symmetric(A), k; debug=(ArpackDebug()))
  @test vals ≈ sort(sort(Avals, by=abs,rev=true)[1:k])
  vals, vecs = eigs(Symmetric(A), k; ritzvec=false, which=:LA)
  @test vals ≈ sort(sort(Avals, rev=true)[1:k])

  vals, vecs = eigs(Float32, Symmetric(A), k; ritzvec=false, which=:LA)
  @test eltype(vals) == Float32
  @test vals ≈ Float32.(sort(sort(Avals, rev=true)[1:k]))

  vals, vecs = eigs(Float32, Float64, Symmetric(A), k; ritzvec=false, which=:LA)
  @test eltype(vals) == Float64
  @test eltype(vecs) == Float32
  @test vals ≈ Float32.(sort(sort(Avals, rev=true)[1:k])) # not too much accuracy

  fop = ArpackSimpleFunctionOp((y,x) -> mul!(y,A,x), n)
  vals, vecs = symeigs(fop, k)
  @test vals ≈ sort(sort(Avals, by=abs,rev=true)[1:k])
end 



@testset "svd interface" begin
  @testset "all ones" begin 
    A = ones(5,5)
    U, s, V = svds(A, 1)
    check_svd(A, U, s, V; tol=10) # not quite as accurate on this one...
    @test s ≈ [5.0]
  end 
  
  @testset "singular case" begin 
    A = ones(10,4)
    @test_throws GenericArpack.ArpackException svds(A, 1; which=:SA, ncv=2)
    
    U,s,V = svds(A, 2; which=:SA, ritzvec=false )
    @test s ≈ [0.0; 0.0] atol=eps(Float64)
    @test size(U,2) == 0
    @test size(V,2) == 0
  end 

  @testset "mytestmat(10,8)" begin 
    A = mytestmat(10,8)
    U, s, V = svds(A, 2; which=:BE)
    check_svd(A, U, s, V; tol=25) # not quite as accurate on this one...
    @test s ≈ [0.20022491452411176, 2.424417285164735]
  end 
end

@testset "generalized interface with tridiag" begin
  n = 100
  A = Tridiagonal(-ones(n-1),2*ones(n),-ones(n-1)).*(n+1)
  #B = Tridiagonal(ones(n-1),4*ones(n),ones(n-1)).*(1/(6*(n+1)))
  B = SymTridiagonal(4*ones(n),ones(n-1)).*(1/(6*(n+1)))
  ABvals,ABvecs = eigen(Matrix(A),Matrix(B))
  k = 4 
  #vals, vecs = eigs(Symmetric(A), k; debug=set_debug_high!(ArpackDebug()))
  vals, vecs = eigs(Symmetric(A), Symmetric(B), k; debug=(ArpackDebug()))
  @test vals ≈ sort(sort(ABvals, by=abs,rev=true)[1:k])

  # needs more iterations...
  vals, vecs = eigs(Symmetric(A), Symmetric(B), k; debug=(ArpackDebug()), which=:SA, maxiter=400, ncv=10)
  @test vals ≈ sort(sort(ABvals)[1:k])
end

@testset "generalized interface" begin 
  using SparseArrays, LinearAlgebra, StableRNGs
  n = 15 
  rng = StableRNG(1234)
  A = sprandn(rng, n, n, 0.25) |> A -> A + A' 
  B = sprandn(rng, n, n, 0.25) |> A -> A*A' 
  ABvals,ABvecs = eigen(Matrix(A),Matrix(B))
  #@show A.colptr, A.rowval, A.nzval

  k = 3 
  #vals, vecs = eigs(Symmetric(A), k; debug=set_debug_high!(ArpackDebug()))
  vals, vecs = eigs(Symmetric(A), Symmetric(B), k; debug=(ArpackDebug()))
  @test vals ≈ sort(sort(ABvals, by=abs,rev=true)[1:k])
  vals, vecs = eigs(Symmetric(A), Symmetric(B), k; ritzvec=false, which=:LA)
  @test vals ≈ sort(sort(ABvals, rev=true)[1:k])

  vals, vecs = eigs(Float32, Symmetric(A), Symmetric(B), k; ritzvec=false, which=:SA)
  @test eltype(vals) == Float32
  @test vals ≈ Float32.(sort(sort(ABvals)[1:k], rev=true))

  #vals, vecs = eigs(Float64, Symmetric(A), k; which=:SM)
  #@test vals ≈ sort(sort(Avals, by=abs)[1:k])
  
  #vals, vecs = eigs(Symmetric(A), k; which=:LA)
end 

@testset "generalized interface with tridiag" begin
  n = 100
  A = Tridiagonal(-ones(n-1),2*ones(n),-ones(n-1)).*(n+1)
  #B = Tridiagonal(ones(n-1),4*ones(n),ones(n-1)).*(1/(6*(n+1)))
  B = SymTridiagonal(4*ones(n),ones(n-1)).*(1/(6*(n+1)))
  ABvals,ABvecs = eigen(Matrix(A),Matrix(B))
  k = 4 
  #vals, vecs = eigs(Symmetric(A), k; debug=set_debug_high!(ArpackDebug()))
  vals, vecs = eigs(Symmetric(A), Symmetric(B), k; debug=(ArpackDebug()))
  @test vals ≈ sort(sort(ABvals, by=abs,rev=true)[1:k])
  vals, vecs = eigs(Symmetric(A), Symmetric(B), k; ritzvec=false, which=:LA)
  @test vals ≈ sort(sort(ABvals, rev=true)[1:k])

  vals, vecs = eigs(Float32, Symmetric(A), Symmetric(B), k; ritzvec=false, which=:SA)
  @test eltype(vals) == Float32
  @test_broken vals ≈ Float32.(sort(sort(ABvals, rev=true)[1:k]))

  vals, vecs = symeigs(Float32, Symmetric(A), Symmetric(B), k; ritzvec=false, which=:LM)
end

@testset "arpack svd example prob" begin
  m = 500 
  n = 100 
  function arpack_svd_example_av!(w,x)
    # computes w = A*x 
    h = 1/(m+1)
    k = 1/(n+1)
    for i=1:m 
      w[i] = 0
    end 
    t = 0.0
    for j=1:n
      t += k
      s = 0.0 
      for i=1:j
        s += h
        w[i] = w[i] + k*s*(t-1)*x[j]
      end
      for i=j+1:m 
        s += h
        w[i] = w[i] + k*t*(s-1)*x[j]
      end
    end
  end
  function arpack_svd_example_atv!(y,w)
    # computes y = A*w
    h = 1/(m+1)
    k = 1/(n+1)
    for i=1:n
      y[i] = 0
    end 
    t = 0.0
    for j=1:n
      t += k 
      s = 0.0
      for i=1:j
        s += h
        y[j] = y[j] + k*s*(t-1)*w[i]
      end 
      for i=j+1:m
        s += h
        y[j] = y[j] + k*t*(s-1)*w[i]
      end 
    end 
  end
  op = ArpackNormalFunctionOp(arpack_svd_example_av!, arpack_svd_example_atv!, m, n)
  einfo = symeigs(op, 4)
  # these are the results from arpack...
  info = svds(op, 4; ncv=10) 
  check_svd(arpack_svd_example_av!, info...)
  @test info.S ≈ sort([4.1012320445852E-02, 6.0488061100249E-02, 1.1784357891005E-01, 5.5723400180223E-01], rev=true)
  info = svds(Float32, op, 4; ncv=10)
  @test info.S ≈ sort([4.10123E-02, 6.04880E-02, 1.17844E-01, 5.57234E-01], rev=true)
  check_svd(arpack_svd_example_av!, info...)
end