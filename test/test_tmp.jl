@testset "tmp" begin
    ss = SimpleSampleSpace(5:25, 75)
    using Gurobi
    solver = GurobiSolver(
        IntFeasTol = 1e-9,
        MIPGapAbs = 1e-3,
        MIPGap = 1e-4,
        Heuristics = .25,
        NumericFocus = 3,
        MIPFocus = 1,
        TimeLimit = 900,
        OutputFlag = 0
    )
    # params = MinimalMinimalExpectedSampleSize(
    #     ss, .2, .4, .05, .2, .4
    # )
    # design, res = getoptimaldesign(params, solver, VERBOSE = 0)
    # println(convert(DataFrames.DataFrame, design))
    # params = SimpleMinimalExpectedSampleSize(
    #     ss, .2, .4, .05, .2, .4
    # )
    # design, res = getoptimaldesign(params, solver, VERBOSE = 0)
    # println(convert(DataFrames.DataFrame, design))

    # params = SimpleMinimalExpectedSampleSize(
    #     ss, .2, .4, .05, .2, .4, GROUPSEQUENTIAL = true
    # )
    # design, res = getoptimaldesign(params, solver, VERBOSE = 1)
    # println(convert(DataFrames.DataFrame, design))
    # println(stoppingforfutility(design, .2))

    params = SimpleMinimalExpectedSampleSize(
        ss, .2, .4, .05, .2, .4, GROUPSEQUENTIAL = true, minstoppingforfutility = .75
    )
    design, res = getoptimaldesign(params, solver, VERBOSE = 1)
    println(convert(DataFrames.DataFrame, design))
    println(stoppingforfutility(design, .2))
    println(score(design))

    params = SimpleMinimalExpectedSampleSize(
        ss, .2, .4, .05, .2, .4, minstoppingforfutility = .75
    )
    design2, res = getoptimaldesign(params, solver, VERBOSE = 1)
    println(convert(DataFrames.DataFrame, design2))
    println(stoppingforfutility(design2, .2))
    println(score(design2))
end
