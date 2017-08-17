"""
    CPInterval <: ConfidenceInterval

    CPInterval{T<:Real}(
        design::Design;
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
julia> ci = CPInterval(est, confidence = .9)
```
"""
struct CPInterval{TD<:Design,TR<:Real} <: ConfidenceInterval

  design::TD
  confidence::TR

  function CPInterval{TD,TR}(
    design::TD, confidence::TR
  ) where {TD<:Design,TR<:Real}

    @checkprob confidence
    new(design, confidence)

  end # inner constructor

end # CPInterval


function CPInterval(
  design::TD, confidence::TR = .9
) where {TD<:Design,TR<:Real}

  return CPInterval{TD,TR}(design::TD, confidence::TR)
  
end # CPInterval


function limits(ci::CPInterval, x1::T, x2::T) where {T<:Integer}

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
