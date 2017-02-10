@testset "MLE" begin
    ss = SimpleSampleSpace(5:30, 75)
    using Gurobi
    solver = GurobiSolver(
        IntFeasTol = 1e-9,
        MIPGapAbs = 1e-3,
        MIPGap = 1e-4,
        Heuristics = .25,
        NumericFocus = 3,
        Threads = 2,
        MIPFocus = 1,
        TimeLimit = 900
    )

    params = SimpleMinimalExpectedSampleSize(
        ss, .2, .4, .05, .2, .4
    )
    design = getoptimaldesign(15, params, solver)
    mle    = MaximumLikelihoodEstimator(design)
    ube    = RaoBlackwellizedEstimator(design)
    cpe    = CompatibleEstimator(design, solver)
    supp   = support(design)
    @test all(
        estimate.(mle, supp[:, 1], supp[:, 2]) .==
            (supp[:, 1] + supp[:, 2])./samplesize.(design, supp[:, 1])
    )
    @test_approx_eq_eps(maximum(abs(bias.(ube, linspace(0, 1)))), 0.0, 1e-10)
    prinln("how do i test this?")
    println(bias.(cpe, linspace(0, 1, 10)))
    println(rmse.(cpe, linspace(0, 1, 10)))
end
