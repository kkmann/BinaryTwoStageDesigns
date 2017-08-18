@testset "'MESS' ............................... " begin

  ss = SampleSpace(5:20, 100)

  params = MESS(
      ss, .2, .4, alpha = .05, beta = .2, pess = .4
  )

  @test null(params) == .2

  @test mtoer(params) == .05

  @test samplespace(params) == ss

  @test maxsamplesize(params) == 100

  @test minconditionalpower(params) == 0.0

  @test typeof(params) == MESS{typeof(1.0)}

  params = MESS(
      ss, .2, .4, minconditionalpower = .7
  )
  @test minconditionalpower(params) == 0.7

  design = optimaldesign(10, params, solver)

  @test score(design) < 26.5

end
