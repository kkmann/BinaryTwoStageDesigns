@testset "Estimation ........................... " begin

  function simonsDesign(r1, n1, r, n)
    nvec = [[n1 for x1 in 0:r1]; [n for x1 in (r1 + 1):n1]]
    cvec = [[Inf for x1 in 0:r1]; [r for x1 in (r1 + 1):n1]]
    return Design(nvec, cvec)
  end

  # Simon's designs for beta = .2, alpha = .05, p1 = p0 +0.2
  # cf. Simon, R "Optimal Two-Stage Designs for Phase II Clinical Trials",
  # Controlled Clinical Trials 10:1-10 (1989). (p. 4)
  sd = simonsDesign( 7, 16, 23, 46) # p0 = 0.4

  mle    = MLEstimator(sd)
  ube    = RBEstimator(sd)
  supp   = support(sd)
  @test all(
    estimate.(mle, supp[:, 1], supp[:, 2]) .==
      (supp[:, 1] + supp[:, 2])./samplesize.(sd, supp[:, 1])
  )
  @test maximum(abs.(bias.(ube, linspace(0, 1)))) ≈ 0.0 atol=1e-10

  ss = SampleSpace(5:20, 50)

  params = MESS(
    ss, .2, .4, alpha = .05, beta = .2, pess = .4
  )

  design = optimaldesign(10, params, solver)

  cpe    = OCEstimator(design, qsolver)

  @test maximum(rmse.(cpe, linspace(0, 1, 10))) <= 0.12

  @test length(incompatibleoutcomes(cpe)) == 0

end
