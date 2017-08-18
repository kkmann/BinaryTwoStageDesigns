@testset "'Design' ............................. " begin

  srand(1234)

  function simonsDesign(r1, n1, r, n)
    nvec = [[n1 for x1 in 0:r1]; [n for x1 in (r1 + 1):n1]]
    cvec = [[Inf for x1 in 0:r1]; [r for x1 in (r1 + 1):n1]]
    return Design(nvec, cvec)
  end

  # Simon's designs for beta = .2, alpha = .05, p1 = p0 +0.2
  # cf. Simon, R "Optimal Two-Stage Designs for Phase II Clinical Trials",
  # Controlled Clinical Trials 10:1-10 (1989). (p. 4)
  p0 = collect(linspace(.1, .7, 7))
  sd = [
    simonsDesign( 1, 10,  5, 29), # p0 = 0.1
    simonsDesign( 3, 13, 12, 43), # p0 = 0.2
    simonsDesign( 5, 15, 18, 46), # p0 = 0.3
    simonsDesign( 7, 16, 23, 46), # p0 = 0.4
    simonsDesign( 8, 15, 26, 43), # p0 = 0.5
    simonsDesign( 7, 11, 30, 43), # p0 = 0.6
    simonsDesign( 4,  6, 22, 27)  # p0 = 0.7
  ]


  
  @test typeof(sd[1]) == Design{typeof(1), typeof(Inf), NoParameters}

  @test DataFrames.DataFrame(sd[1]) == DataFrames.DataFrame(
        x1 = 0:interimsamplesize(sd[1]),
        n  = samplesize(sd[1]),
        c  = criticalvalue(sd[1])
    )

  @test interimsamplesize(sd[1]) == 10

  @test parameters(sd[1]) == NoParameters()

  @test samplesize(sd[1]) == [[10 for x1 in 0:1]; [29 for x1 in (1 + 1):10]]
  @test samplesize(sd[1], 1) == 10
  @test samplesize(sd[1], 2) == 29

  @test criticalvalue(sd[1]) == [[Inf for x1 in 0:1]; [5 for x1 in (1 + 1):10]]
  @test criticalvalue(sd[1], 1) == Inf
  @test criticalvalue(sd[1], 2) == 5

  supp_ = support(sd[1])
  for p in linspace(0, 1, 11)
    
    @test pdf(sd[1], 1, 2, p) == 0.0
    pdf_ = pdf.(sd[1], support(sd[1])[:,1], support(sd[1])[:,2], p)
    @test all(pdf_ .>= 0.0)
    @test isapprox(sum(pdf_), 1.0, atol = 0.00001)

  end

  @test power(sd[1], p0[1]) <= .05
  @test power(sd[1], p0[1] + .2) >= .8

  prior = p -> 1
  @test expectedpower(sd[1], prior, mcrv = p0[1] + .1) < power(sd[1], .8)
  @test expectedpower(sd[1], prior, mcrv = p0[1] + .1) > power(sd[1], .2)
  @test expectedpower(sd[1], 1, prior, mcrv = p0[1] + .1) < expectedpower(sd[1], 2, prior, mcrv = p0[1] + .1)

  @test isapprox(
    stoppingforfutility(sd[1], .2), pdf(sd[1], 0, 0, .2) + pdf(sd[1], 1, 0, .2), atol = .00001
  )

  @test !test(sd[1], 0, 0)
  @test test(sd[1], 5, 11)

  sim_ = simulate(sd[1], p0[1], 1000)
  @test mean(sim_[:rejectedH0]) < .05
  @test all(sim_[:n] .<= 29)
  @test all(sim_[:n] .>= 0)

  @test isapprox(quadgk(jeffreysprior(sd[1]), 0, 1)[1], 1.0, atol = 0.0001)

  save("test.jls", sd[1])
  @test true
  rm("test.jls")
  writecsv("test.csv", sd[1]; label = "test")
  @test true
  rm("test.csv")

end
