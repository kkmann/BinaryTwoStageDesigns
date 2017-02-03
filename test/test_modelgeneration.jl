@testset "model generation" begin

    ss = SimpleSampleSpace(5:30, 125)
    using Gurobi

    params = SimpleMinimalExpectedSampleSize(
        ss, .2, .4, .05, .2, .4, Unimodal, StoppingForEfficacy
    )
    m, y = BinaryTwoStageDesigns._createProblem(10, params)
    getoptimaldesign(20, params, GurobiSolver())

end
