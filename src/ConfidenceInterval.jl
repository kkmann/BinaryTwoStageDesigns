abstract BinaryTwoStageDesignConfidenceInterval

# these two enable array broadcasting of all methods!
Base.size(::BinaryTwoStageDesignConfidenceInterval) = ()
Base.getindex(ci::BinaryTwoStageDesignConfidenceInterval, i) = ci

function limits{T<:Integer}(ci::BinaryTwoStageDesignConfidenceInterval, x1::T, x2::T)
    error("confidence interval is not implemented!")
end

function getConfidence(ci::BinaryTwoStageDesignConfidenceInterval)
    try
        return ci.confidence
    catch
        error("not implemented")
    end
end
