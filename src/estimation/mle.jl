type MaximumLikelihoodEstimator <: BinaryTwoStageDesignEstimator
    design::AbstractBinaryTwoStageDesign

    function MaximumLikelihoodEstimator(
        design::AbstractBinaryTwoStageDesign
    )
        new(design)
    end
end

function estimate{T<:Integer}(estimator::MaximumLikelihoodEstimator, x1::T, x2::T)::Float64
    supp = _support(estimator.design)
    if !_isInSupport(supp, x1, x2)
        return NaN
    else
        return (x1 + x2)/getSampleSize(estimator.design, x1)
    end
end
