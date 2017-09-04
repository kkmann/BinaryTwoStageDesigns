@testset "'LiuScore' ........................... " begin

  ss = SampleSpace(5:20, 60)
  prior = p -> Distributions.pdf(Distributions.Beta(2 , 5), p)
  params = LiuScore(ss, .2, .3, prior, .05, .2, 2.0, 0.2)
  design = optimaldesign(10, params, solver)
  @test power(design, .2) <= .05

end
