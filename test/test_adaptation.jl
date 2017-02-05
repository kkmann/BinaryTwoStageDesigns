@testset "design adaptation" begin

    ss = SimpleSampleSpace(5:30, 75)
    using Gurobi
    solver = GurobiSolver(
        IntFeasTol = 1e-9,
        MIPGapAbs = 1e-3,
        MIPGap = 1e-4,
        Heuristics = .25,
        NumericFocus = 3,
        # Threads = 2,
        MIPFocus = 1,
        TimeLimit = 900
    )

    params = SimpleMinimalExpectedSampleSize(
        ss, .2, .4, .05, .2, .4
    )
    design = getoptimaldesign(15, params, solver)
    println(convert(DataFrames.DataFrame, design))

    outcome = DataArrays.@data(
        [true; true; true; true; false;
         false; false; false; false; false;
         false; false; false; NA]
    )
    stage = repmat([1], 14)
    adapteddesign = adapt(design, outcome, stage, solver)
    println(convert(DataFrames.DataFrame, adapteddesign))

end
