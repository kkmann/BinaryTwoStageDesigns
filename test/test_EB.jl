@testset "'EB' ................................. " begin

  ss = SampleSpace(5:20, 60)
  prior = p -> Distributions.pdf(Distributions.Beta(2 , 5), p)
  params = EB(ss, .2, .3, prior, 16*10^6, 50*10^3)
  design = optimaldesign(10, params, solver)
  @test power(design, .2) <= .05

  params = EB(ss, .2, .3, prior, 16*10^6, 50*10^3, a = 21, b = 9)
  design = optimaldesign(10, params, solver)
  @test power(design, .2) <= .05

  params = EB(ss, .2, .3, prior, 16*10^6, 50*10^3, a = Inf, b = Inf, targetpower = 0.8)
  design = optimaldesign(10, params, solver)
  @test power(design, .2) <= .05

end
