@testset "tmp" begin
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
    params = MinimalMinimalExpectedSampleSize(
        ss, .2, .4, .05, .2, .4
    )
    design, res = getoptimaldesign(params, solver, VERBOSE = 0)
    println(convert(DataFrames.DataFrame, design))
    params = SimpleMinimalExpectedSampleSize(
        ss, .2, .4, .05, .2, .4
    )
    design, res = getoptimaldesign(params, solver, VERBOSE = 0)
    println(convert(DataFrames.DataFrame, design))
end
