@testset "'SampleSpace' ........................ " begin

  ss = SampleSpace(collect(5:10), 25)

  ss = SampleSpace(5:10, 25)

  @test interimsamplesizerange(ss) == collect(5:10)

  @test maxsamplesize(ss) == 25

  @test all(possible.(5:10, ss))

  @test all(.!possible.([-3, 0, 4, 11, 20], ss))

  @test isgroupsequential(ss) == false

  ss = SampleSpace(5:10, 25, n2min = 10)
  ss = SampleSpace(5:10, 25, maxnfact = 5)
  ss = SampleSpace(5:10, 25, nmincont = 10)
  ss = SampleSpace(5:10, 25, maxvariables = 10000)

end
