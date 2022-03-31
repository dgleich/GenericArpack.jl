@testset "dsaupd with idonow" begin 
  using LinearAlgebra
  op = ArpackInJulia.ArpackSimpleOp(Diagonal(1.0:30))
  prob1 = _allocate_symproblem(op, 12)
  prob2 = deepcopy(prob1)

  niter1 = _eigrun!(prob1, 6)
  niter2 = _eigrun!(prob2, 6; idonow=true)

  @test niter2 == 1
  @test niter1 > niter2 # needs at least one matvec...
  _compare_probs(prob1, prob2) # run tests on each field
  @test prob2.op == prob1.op
end 