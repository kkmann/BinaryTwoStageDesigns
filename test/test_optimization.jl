@testset "'optimaldesign' ...................... " begin

  ss = SampleSpace(15:17, 60)
  params = MESS(ss, .2, .4)
  design, res = optimaldesign(params, solver, EARLYTERMINATION = true, VERBOSE = 0)

  @test all(score(design) .<= res["scores"])

end
