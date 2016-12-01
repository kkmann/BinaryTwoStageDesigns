type ClopperPearsonConfidenceInterval <: BinaryTwoStageDesignConfidenceInterval
    estimator::BinaryTwoStageDesignEstimator
    confidence::Float64

    function ClopperPearsonConfidenceInterval{T<:Real}(
        estimator::BinaryTwoStageDesignEstimator;
        confidence::T = .9
    )
        @assert 0 < confidence
        @assert confidence < 1
        new(estimator, confidence)
    end
end

function limits{T<:Integer}(ci::ClopperPearsonConfidenceInterval, x1::T, x2::T; k::Integer = 1001)
    @assert _isInSupport(_support(ci.estimator.design), x1, x2)
    n1     = getInterimSampleSize(ci.estimator.design)
    nmax   = maximum(getInterimSampleSize(ci.estimator.design))
    grid   = collect(linspace(0, 1, k))
    pvals  = map(p0 -> p(ci.estimator, x1, x2, p0), grid)
    indlow = findfirst(pvals .> (1 - ci.confidence)/2.0)
    indup  = findlast(pvals .< 1 - (1 - ci.confidence)/2.0)
    limits = [grid[indlow]; grid[indup]]
    if indlow == 1
        limits[1] = 0.0
    end
    if indup == k
        limits[2] = 1.0
    end
    if x1 + x2 == 0
        limits[1] = 0.0
    end
    if x1 + x2 == getSampleSize(ci.estimator.design, x1)
        limits[2] = 1.0
    end
    return limits
end
