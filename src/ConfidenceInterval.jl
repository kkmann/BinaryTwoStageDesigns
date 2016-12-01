abstract BinaryTwoStageDesignConfidenceInterval

function limits{T<:Integer}(ci::BinaryTwoStageDesignConfidenceInterval, x1::T, x2::T)
    error("confidence interval is not implemented!")
end

function getConfidence(ci::BinaryTwoStageDesignConfidenceInterval)
    try
        return estimator.confidence
    catch
        error("not implemented")
    end
end
