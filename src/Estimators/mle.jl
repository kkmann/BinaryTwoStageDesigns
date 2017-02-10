type MaximumLikelihoodEstimator <: BinaryTwoStageDesignEstimator
    design::AbstractBinaryTwoStageDesign

    function MaximumLikelihoodEstimator(
        design::AbstractBinaryTwoStageDesign
    )
        new(design)
    end
end

function estimate{T<:Integer}(estimator::MaximumLikelihoodEstimator, x1::T, x2::T)
    checkx1x2(x1, x2, design(estimator))
    return (x1 + x2)/samplesize(design(estimator), x1)
end
