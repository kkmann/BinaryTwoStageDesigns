@testset "model generation" begin

    ss = SimpleSampleSpace(5:30, 125)
    using Gurobi
    solver = GurobiSolver(
        IntFeasTol = 1e-9,
        MIPGapAbs = 1e-3,
        MIPGap = 1e-5,
        Heuristics = .2,
        NumericFocus = 3
    )
    params = SimpleMinimalExpectedSampleSize(
        ss, .2, .4, .05, .2, .4, Unimodal, StoppingForEfficacy
    )
    m, y = BinaryTwoStageDesigns._createProblem(10, params)
    getoptimaldesign(20, params, solver)

end
