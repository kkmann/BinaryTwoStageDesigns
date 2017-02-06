@testset "check simons designs" begin

p0 = collect(linspace(.1, .7, 7))
p1 = p0 + .2
r1 = [1 3 5 7 8 7 4]
n1 = [10 13 15 16 15 11 6]
r  = [5 12 18 23 26 30 22]
n  = [29 43 46 46 43 43 27]

import Gurobi

solver = Gurobi.GurobiSolver(
    IntFeasTol   = 1e-9,
    MIPGapAbs    = 1e-3,
    MIPGap       = 1e-4,
    Heuristics   = .25,
    NumericFocus = 3,
    Threads      = 2,
    MIPFocus     = 1,
    TimeLimit    = 900,
    OutputFlag   = 0
)
ss = SimpleSampleSpace(3:20, 70)

function simonsdesign(r1, n1, r, n)
    nvec = [[n1 for x1 in 0:r1]; [n for x1 in (r1 + 1):n1]]
    cvec = [[Inf for x1 in 0:r1]; [r for x1 in (r1 + 1):n1]]
    return BinaryTwoStageDesign(nvec, cvec)
end

    @testset "compare with optimization results" for i in 1:length(p0)
        sd = simonsdesign( r1[i], n1[i], r[i], n[i])
        params = SimpleMinimalExpectedSampleSize(
            ss, p0[i], p1[i], .05, .2, p0[i], GROUPSEQUENTIAL = true, STOPPINGFOREFFICACY = false
        )
        d, other = getoptimaldesign(params, solver)
        @test samplesize(d) == samplesize(sd)
        @test criticalvalue(d) == criticalvalue(sd)
    end
end
