@testset "SampleSize" begin

    function simonsDesign(r1, n1, r, n)
        nvec = [[n1 for x1 in 0:r1]; [n for x1 in (r1 + 1):n1]]
        cvec = [[Inf for x1 in 0:r1]; [r for x1 in (r1 + 1):n1]]
        return BinaryTwoStageDesign(nvec, cvec)
    end

    # Simon's designs for beta = .2, alpha = .05, p1 = p0 +0.2
    # cf. Simon, R "Optimal Two-Stage Designs for Phase II Clinical Trials",
    # Controlled Clinical Trials 10:1-10 (1989). (p. 4)
    p0 = collect(linspace(.1, .7, 7))
    sd = [
        simonsDesign( 1, 10,  5, 29), # p0 = 0.1
        simonsDesign( 3, 13, 12, 43), # p0 = 0.2
        simonsDesign( 5, 15, 18, 46), # p0 = 0.3
        simonsDesign( 7, 16, 23, 46), # p0 = 0.4
        simonsDesign( 8, 15, 26, 43), # p0 = 0.5
        simonsDesign( 7, 11, 30, 43), # p0 = 0.6
        simonsDesign( 4,  6, 22, 27)  # p0 = 0.7
    ]


    println("compare expected sample size with publication ... ")
    en = [15.0, 20.6, 23.6, 24.5, 23.5, 20.5, 14.8]
    for i in 1:length(p0)
        ss = SampleSize(sd[i], p0[i])
        @test abs(mean(ss) - en[i]) < .1
        # should be sufficient to ensure correctnes,
        # check whether all routines can be called
        nr = rand(ss, convert(Int64, 1e4))
        @test all(minimum(ss) .<= nr .<= maximum(ss))
        pdf(ss, 25)
        quantile(ss, .5)
    end

end
