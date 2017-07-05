"""
    RaoBlackwellizedEstimator <: BinaryTwoStageDesignEstimator

    RaoBlackwellizedEstimator(design::BinaryTwoStageDesign)

Unbiased estimator for response rate `p` see also:

Kunzmann K, Kieser M. Point estimation and p‐values in phase II adaptive two‐stage designs with a binary endpoint. Statistics in medicine. 2017 Mar 15;36(6):971-84.

# Examples
```julia-repl
julia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)
julia> interimsamplesize(ss)
julia> design = getoptimaldesign(15, params, solver = Gurobi.GurobiSolver())
julia> est = RaoBlackwellizedEstimator(design)
```
"""
type RaoBlackwellizedEstimator <: BinaryTwoStageDesignEstimator
    design::BinaryTwoStageDesign

    function RaoBlackwellizedEstimator(
        design::BinaryTwoStageDesign
    )
        new(design)
    end
end

function estimate{T<:Integer}(estimator::RaoBlackwellizedEstimator, x1::T, x2::T)
    checkx1x2(x1, x2, design(estimator))
    n1          = interimsamplesize(design(estimator))
    nominator   = BigInt(0.0) # we need big ints, as binomials get HUGE!
    denominator = BigInt(0.0)
    x = x1 + x2
    n = samplesize(design(estimator), x1)
    for xx1 in 0:(min(x, n1))
        if samplesize(design(estimator), xx1) == n
            nominator   += BigInt(binomial(n1 - 1, xx1 - 1))*BigInt(binomial(n - n1, x - xx1))
            denominator += BigInt(binomial(n1, xx1))*BigInt(binomial(n - n1, x - xx1))
        end
    end
    return convert(Float64, nominator/denominator)
end
