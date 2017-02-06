@testset "testing simple minimal expected sample size type" begin

    ss = SimpleSampleSpace(5:20, 100)

    params = SimpleMinimalExpectedSampleSize(
        ss, .2, .4, .05, .2, .4
    )
    @test null(params) == .2
    @test alpha(params) == .05
    @test samplespace(params) == ss
    @test maxsamplesize(params) == 100
    @test allowsstoppingforefficacy(params) == true
    @test isgroupsequential(params) == false
    @test minconditionalpower(params) == 0.0
    @test typeof(params) == SimpleMinimalExpectedSampleSize{
        SimpleSampleSpace{typeof(100)},
        NotGroupSequential,
        AllowStoppingForEfficacy,
        NoMonotoneConditionalPower
    }

    params = SimpleMinimalExpectedSampleSize(
        ss, .2, .4, .05, .2, .4, minconditionalpower = .7, GROUPSEQUENTIAL = true, STOPPINGFOREFFICACY = false
    )
    @test allowsstoppingforefficacy(params) == false
    @test isgroupsequential(params) == true
    @test minconditionalpower(params) == 0.7
    @test typeof(params) == SimpleMinimalExpectedSampleSize{
        SimpleSampleSpace{typeof(100)},
        GroupSequential,
        NoStoppingForEfficacy,
        NoMonotoneConditionalPower
    }
end

# TODO test optimization
