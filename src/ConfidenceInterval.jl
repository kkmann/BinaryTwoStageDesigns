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

function getDesign(ci::BinaryTwoStageDesignConfidenceInterval)
    try
        return ci.design
    catch
        error("not implemented")
    end
end

function coverage{T<:Real}(
    ci::BinaryTwoStageDesignConfidenceInterval,
    p::T
)
    @assert 0 <= p
    @assert p <= 1
    design = getDesign(ci)
    n1     = getInterimSampleSize(design)
    nmax   = maximum(getSampleSize(design))
    # construct array of possible outcomes
    supp   = _support(design)
    nposs  = size(supp, 1)
    res    = 0.0
    for i in 1:nposs
        x1 = supp[i, 1]
        x2 = supp[i, 2]
        lim = limits(ci, x1, x2)
        if (lim[1] <= p) & (lim[2] >= p)
            res += probability(design, x1, x2, p)
        end
    end
    return res
end
