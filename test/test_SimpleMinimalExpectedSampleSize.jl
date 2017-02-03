@testset "SimpleMinimalExpectedSampleSize" begin

    ss = SimpleSampleSpace(5:20, 100)

    params = SimpleMinimalExpectedSampleSize(
        ss, .2, .4, .05, .2, .4, GroupSequential, StoppingForEfficacy
    )
    @test null(params) == .2
    @test alpha(params) == .05
    @test samplespace(params) == ss
    @test maxsamplesize(params) == 100
    @test typeof(params) == SimpleMinimalExpectedSampleSize{SimpleSampleSpace{typeof(100)}, GroupSequential, StoppingForEfficacy}
    @test efficacy(params) == StoppingForEfficacy
    @test regularization(params) == GroupSequential
end

# TODO test optimization
