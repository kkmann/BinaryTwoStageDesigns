@testset "testing simple sample space type" begin

    ss = SimpleSampleSpace(collect(5:10), 25)

    ss = SimpleSampleSpace(5:10, 25)

    @test interimsamplesizerange(ss) == collect(5:10)
    @test maxsamplesize(ss) == 25
    @test all(possible.(5:10, ss))
    @test all(!possible.([-3, 0, 4, 11, 20], ss))

end
