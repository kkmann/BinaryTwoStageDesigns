@testset "testing simple minimal expected sample size type" begin

    ss = SimpleSampleSpace(5:20, 100)

    params = SimpleMinimalExpectedSampleSize(
        ss, .2, .4, .05, .2, .4
    )
    @test null(params) == .2
    @test alpha(params) == .05
    @test samplespace(params) == ss
    @test maxsamplesize(params) == 100
    @test minconditionalpower(params) == 0.0
    @test typeof(params) == SimpleMinimalExpectedSampleSize{SimpleSampleSpace{typeof(100)}}

    params = SimpleMinimalExpectedSampleSize(
        ss, .2, .4, .05, .2, .4, minconditionalpower = .7
    )
    @test minconditionalpower(params) == 0.7
    @test typeof(params) == SimpleMinimalExpectedSampleSize{SimpleSampleSpace{typeof(100)}}
end
