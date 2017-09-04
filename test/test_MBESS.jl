@testset "'MBESS' .............................. " begin

  ss = SampleSpace(5:20, 60)
  prior = p -> Distributions.pdf(Distributions.Beta(2 , 5), p)
  params = MBESS(ss, .2, .3, prior, a = 1, b = 1)
  design = optimaldesign(10, params, solver)

  @test expectedpower(design, prior) > .2

  params = MBESS(ss, .2, .3, prior, a = Inf, b = Inf, targetpower = 0.8)
  design = optimaldesign(10, params, solver)

  z =  z = quadgk(prior, mcrv(parameters(design)), 1)[1]
  tmp = fzero(p -> power(design, p) - .8, 0, 1)
  @test quadgk(prior, tmp, 1)[1]/z > 0.8

end
