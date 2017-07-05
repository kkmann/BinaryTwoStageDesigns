"""
    NaiveClopperPearsonConfidenceInterval <: ConfidenceInterval

    NaiveClopperPearsonConfidenceInterval{T<:Real}(
        design::BinaryTwoStageDesign;
        confidence::T = .9
    )

Naive Clopper-Pearson confidence interval using default ordering.

# Parameters

| Parameter    | Description |
| -----------: | :---------- |
| estimator    | estimator object defining the sample space ordering |
| confidence   | confidence level of the interval |

# Examples
```julia-repl
julia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)
julia> interimsamplesize(ss)
julia> design = getoptimaldesign(15, params, Gurobi.GurobiSolver())
julia> est = MaximumLikelihoodEstimator(design, Gurobi.GurobiSolver())
julia> ci = NaiveClopperPearsonConfidenceInterval(est, confidence = .9)
```
"""
immutable NaiveClopperPearsonConfidenceInterval <: ConfidenceInterval
    design::BinaryTwoStageDesign
    confidence::Float64

    function NaiveClopperPearsonConfidenceInterval{T<:Real}(
        design::BinaryTwoStageDesign;
        confidence::T = .9
    )
        checkp(confidence)
        new(design, confidence)
    end
end #NaiveClopperPearsonConfidenceInterval



function limits{T<:Integer}(ci::NaiveClopperPearsonConfidenceInterval, x1::T, x2::T)
    d              = design(ci)
    ispossible(d, x1, x2) ? nothing : error("(x1, x2) not compatible with underlying design")
    alpha::Float64 = 1 - confidence(ci)
    n1             = interimsamplesize(d)
    nmax::Int64    = maximum(samplesize(d))
    x              = x1 + x2
    n              = samplesize(d, x1)
    if x == 0
        limits = [0.0; (1.0 - (alpha/2)^(1/n))]
    elseif x == n
        limits = [(alpha/2)^(1/n); 1.0]
    else
        limits = [quantile(Distributions.Beta(x, n - x + 1), alpha/2); quantile(Distributions.Beta(x + 1, n - x), 1 - alpha/2)]
    end
    return limits
end
