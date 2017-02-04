@testset "model generation" begin

    ss = SimpleSampleSpace(5:30, 75)
    using Gurobi
    solver = GurobiSolver(
        IntFeasTol = 1e-9,
        MIPGapAbs = 1e-3,
        MIPGap = 1e-4,
        Heuristics = .25,
        NumericFocus = 3,
        Threads = 2,
        MIPFocus = 1
        #OutputFlag = 0 # TODO might not be a good idea...
    )
    params = SimpleMinimalExpectedSampleSize(
        ss, .2, .4, .05, .2, .4, Unimodal, StoppingForEfficacy
    )
    m, y = BinaryTwoStageDesigns._createProblem(10, params)
    design = getoptimaldesign(15, params, solver)
    println(convert(DataFrames.DataFrame, design))
    #getoptimaldesign(params, solver)

end
