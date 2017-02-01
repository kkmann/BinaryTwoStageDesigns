using Base.Test
using BinaryTwoStageDesigns

@testset "BinaryTwoStageDesigns" begin

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
        simonsDesign( 4,  6, 22, 27) # p0 = 0.7
    ]

    @testset "show() and convert()" begin
        println("testing whether show() and convert() methods can be called ... ")
        originalSTDOUT = STDOUT
        redirect_stdout()
        function testShow(design)
            show(design)
            return true
        end
        for i in 1:length(p0)
            @test testShow(sd[i])
        end
        redirect_stdout(originalSTDOUT)
        @test convert(DataFrames.DataFrame, sd[1]) == DataFrames.DataFrame(
            x1 = 0:getInterimSampleSize(sd[1]),
            n  = getSampleSize(sd[1]),
            c  = getRejectionBoundary(sd[1])
        )
    end

    @testset "Probability to Reject" begin
        println("verifying Simon's designs alpha, beta ... ")
        for i in 1:length(p0)
            @test probabilityToReject(sd[i], p0[i]) <= .05
            @test probabilityToReject(sd[i], p0[i] + 0.2) >= .8
        end
    end

    @testset "get methods" begin
        println("testing get methods ... ")
        @test getInterimSampleSize(sd[1]) == 10
        n = [[10 for x1 in 0:1]; [29 for x1 in 2:10]]
        c = [[Inf for x1 in 0:1]; [5 for x1 in 2:10]]
        @test getSampleSize(sd[1]) == n
        @test getRejectionBoundary(sd[1]) == c
        for x1 in 0:getInterimSampleSize(sd[1])
            @test getSampleSize(sd[1], x1) == n[x1 + 1]
            @test getRejectionBoundary(sd[1], x1) == c[x1 + 1]
        end
    end

    @testset "test()" begin
        println("testing test decisions ... ")
        for x1 in 0:getInterimSampleSize(sd[1])
            for x2 in 0:(getSampleSize(sd[1], x1) - getInterimSampleSize(sd[1]))
                if x1 + x2 > getRejectionBoundary(sd[1], x1)
                    @test test(sd[1], x1, x2)
                else
                    @test !test(sd[1], x1, x2)
                end
            end
        end
    end

    @testset "simulate Simon's designs" begin
        println("validate probability to reject by simulation, potential sampling errors! ... ")
        for i in 1:length(p0)
            simres = simulate(sd[i], p0[i], convert(Int, 1e4))
            psim = mean(simres[:rejectedH0])
            pnum = probabilityToReject(sd[i], p0[i])
            @test_approx_eq_eps(psim, pnum, p0[i]/sqrt(1e4)*5) # 5 sigma should be fine most of the time
            ss = SampleSize(sd[i], p0[i])
            munhat = mean(simres[:n])
            mun    = mean(ss)
            @test_approx_eq_eps(munhat, mun, mun*.05)
            varnhat = var(simres[:n])
            varn    = var(ss)
            @test_approx_eq_eps(varnhat, varn, varn*.05)

            simres = simulate(sd[i], p0[i] + .2, convert(Int, 1e4))
            psim = mean(simres[:rejectedH0])
            pnum = probabilityToReject(sd[i], p0[i] + .2)
            @test_approx_eq_eps(psim, pnum, (p0[i] + .2)/sqrt(1e4)*5) # 5 sigma should be fine most of the time
            ss = SampleSize(sd[i], p0[i] + .2)
            munhat = mean(simres[:n])
            mun    = mean(ss)
            @test_approx_eq_eps(munhat, mun, mun*.05)
            varnhat = var(simres[:n])
            varn    = var(ss)
            @test_approx_eq_eps(varnhat, varn, varn*.05)
        end

    end

    println("done.")

end
