"""
    MaximumLikelihoodEstimator

    MaximumLikelihoodEstimator(design::BinaryTwoStageDesign)

Simple maximum likelihood estimator for response rate `p`.

# Examples
```julia-repl
julia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)
julia> interimsamplesize(ss)
julia> design = getoptimaldesign(15, params, solver = Gurobi.GurobiSolver())
julia> est = MaximumLikelihoodEstimator(design)
```
"""
type MaximumLikelihoodEstimator <: BinaryTwoStageDesignEstimator
    design::BinaryTwoStageDesign

    function MaximumLikelihoodEstimator(
        design::BinaryTwoStageDesign
    )
        new(design)
    end
end

function estimate{T<:Integer}(estimator::MaximumLikelihoodEstimator, x1::T, x2::T)
    checkx1x2(x1, x2, design(estimator))
    return (x1 + x2)/samplesize(design(estimator), x1)
end
