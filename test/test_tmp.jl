@testset "tmp" begin
    ss = SimpleSampleSpace(5:15, 75)
    using Gurobi
    solver = GurobiSolver(
        IntFeasTol = 1e-9,
        MIPGapAbs = 1e-3,
        MIPGap = 1e-4,
        Heuristics = .25,
        NumericFocus = 3,
        MIPFocus = 1,
        TimeLimit = 60,
        OutputFlag = 1
    )
    # params = SimpleMinimalExpectedSampleSize(
    #     ss, .2, .4, .05, .2, .4; minconditionalpower = .7
    # )
    # design, res = getoptimaldesign(params, solver, VERBOSE = 1)
    # println(convert(DataFrames.DataFrame, design))
    # println(score(design))
    #
    # params = SimpleMinimalExpectedSampleSize(
    #     ss, .2, .4, .05, .2, .4, GROUPSEQUENTIAL = true, minconditionalpower = .7
    # )
    # design, res = getoptimaldesign(params, solver, VERBOSE = 1)
    # println(convert(DataFrames.DataFrame, design))
    # println(score(design))

    # params = SimpleMinimalExpectedSampleSize(
    #     ss, .2, .4, .05, .2, .4, GROUPSEQUENTIAL = true, minstoppingforfutility = .75
    # )
    # design, res = getoptimaldesign(params, solver, VERBOSE = 1)
    # println(convert(DataFrames.DataFrame, design))
    # println(stoppingforfutility(design, .2))
    # println(score(design))
    #
    # params = SimpleMinimalExpectedSampleSize(
    #     ss, .2, .4, .05, .2, .4, minstoppingforfutility = .75
    # )
    # design2, res = getoptimaldesign(params, solver, VERBOSE = 1)
    # println(convert(DataFrames.DataFrame, design2))
    # println(stoppingforfutility(design2, .2))
    # println(score(design2))
    function prior(p)
        if (p > .3) & (p < .5)
            return 5.0
        else
            return 0.0
        end
    end
    params = LiuScore(
        ss, .2, prior, .05, .2, .5, .5
    )
    design = getoptimaldesign(15, params, solver)
    # design, res = getoptimaldesign(params, solver, VERBOSE = 1)
    println(convert(DataFrames.DataFrame, design))
    println(score(design))
    println(score.(design, params, linspace(.2, .6, 15)))
    println(rup.(design, params, linspace(.2, .6, 15)))
    println(ros.(design, params, linspace(.2, .6, 15)))
end
