@testset "Confidence intervals ................. " begin
  
  ss = SampleSpace(5:20, 50)

  params = MESS(
    ss, .2, .4, alpha = .05, beta = .2, pess = .4
  )

  design = optimaldesign(10, params, solver)
  supp   = support(design)

  cpe    = OCEstimator(design, qsolver)

  ci     = CPInterval(design, confidence = .9)
  @test minimum(coverage.(ci, linspace(0, 1))) > .85 

  lim  = hcat(limits.(ci, supp[:, 1], supp[:, 2])...)'
  @test all((0 .<= lim[:, 1]) .& (1 .>= lim[:, 2]))

  ci     = ECPInterval(cpe, confidence = .9)
  @test minimum(coverage.(ci, linspace(0, 1))) >= .9
  lim  = hcat(limits.(ci, supp[:, 1], supp[:, 2])...)'
  @test all((0 .<= lim[:, 1]) .& (1 .>= lim[:, 2]))

  # do not test mmwinterval, not recommended for use

end
