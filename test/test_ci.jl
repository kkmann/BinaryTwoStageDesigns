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
    cpe    = CompatibleEstimator(design, solver)
    ci     = ClopperPearsonConfidenceInterval(cpe, confidence = .9)
    println(coverage.(ci, linspace(0, 1)))
    supp = support(design)
    for i in 1:size(supp, 1)
        x1, x2 = supp[i, :]
        lim = limits(ci, x1, x2)
    end
end
