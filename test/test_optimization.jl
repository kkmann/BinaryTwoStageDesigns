@testset "model generation" begin

    ss = SimpleSampleSpace(5:30, 75)
    using Gurobi
    solver = GurobiSolver(
        TimeLimit = 900
    )

    params = SimpleMinimalExpectedSampleSize(
        ss, .2, .4, .05, .2, .4
    )
    design = getoptimaldesign(15, params, solver)
    println(convert(DataFrames.DataFrame, design))
    println(power.(design, 0:15, .4))

    params = SimpleMinimalExpectedSampleSize(
        ss, .2, .4, .05, .2, .4, MONOTONECONDITIONALPOWER = true
    )
    design = getoptimaldesign(15, params, solver)
    println(convert(DataFrames.DataFrame, design))
    println(power.(design, 0:15, .4))

    params = SimpleMinimalExpectedSampleSize(
        ss, .2, .4, .05, .2, .4, minconditionalpower = .7
    )
    design = getoptimaldesign(15, params, solver)
    println(convert(DataFrames.DataFrame, design))

    ss = SimpleSampleSpace(5:30, 75, GS = true)
    params = SimpleMinimalExpectedSampleSize(
        ss, .2, .4, .05, .2, .4
    )
    design = getoptimaldesign(15, params, solver)
    println(convert(DataFrames.DataFrame, design))

    ss = SimpleSampleSpace(5:30, 75, nmincont = 10)
    params = SimpleMinimalExpectedSampleSize(
        ss, .2, .4, .05, .2, .4,
    )
    design = getoptimaldesign(15, params, solver)
    println(convert(DataFrames.DataFrame, design))
end
