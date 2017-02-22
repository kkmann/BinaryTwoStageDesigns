@testset "tmp" begin
    ss = SimpleSampleSpace(20:35, 75)
    using Gurobi
    solver = GurobiSolver(
        IntFeasTol = 1e-9,
        MIPGapAbs = 1e-3,
        MIPGap = 1e-4,
        Heuristics = .25,
        NumericFocus = 3,
        MIPFocus = 1,
        TimeLimit = 120,
        OutputFlag = 1,
        InfUnbdInfo = 1
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
    #     ss, .2, .4, .05, .2, .4, minstoppingforfutility = .7, STOPPINGFOREFFICACY = false
    # )
    # design2, res = getoptimaldesign(params, solver, VERBOSE = 1)
    # println(convert(DataFrames.DataFrame, design2))
    # println(score(design2))
    function prior(p)
        if (p > .3) & (p < .5)
            return 5.0
        else
            return 0.0
        end
    end
    params = LiuScore(
        ss, .2, prior, .05, .2, .5, .5, minconditionalpower = .7
    )
    # design = getoptimaldesign(15, params, solver)
    design, res = getoptimaldesign(params, solver, VERBOSE = 1)
    println(convert(DataFrames.DataFrame, design))
    # println(score(design))
    # priorpivots = collect(linspace(.2, 1.0, 250 + 2))[2:(250 + 1)] # leave out boundary values!
    # dp          = priorpivots[2] - priorpivots[1]
    # priorvals   = prior.(priorpivots) ./ sum(prior.(priorpivots) .* dp) # normalize to 1
    # println(score.(design, params, priorpivots))
    # println(  rup.(design, params, priorpivots))
    # println(  ros.(design, params, priorpivots))
    # println(sum(score.(design, params, priorpivots) .* priorvals .* dp))
    println(score(design))
    println(power.(design, 0:interimsamplesize(design), .4))
end
