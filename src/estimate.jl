"""
Abstract base type for estimators
"""
abstract BinaryTwoStageDesignEstimator

function getDesign(estimator::BinaryTwoStageDesignEstimator)
    try
        return estimator.design
    catch
        error("not implemented")
    end
end

function estimate{T<:Integer}(estimator::BinaryTwoStageDesignEstimator, x1::T, x2::T)::Float64
    error("estimation is not implemented!")
end

"""
Return natural p-value for given estimator (H0: p <= p0)
"""
function p{T<:Integer}(estimator::BinaryTwoStageDesignEstimator, x1::T, x2::T)::Float64
    return p(estimator, x1, x2, estimator.design.parameters.p0)
end
function p{T1<:Integer, T2<:Real}(estimator::BinaryTwoStageDesignEstimator, x1::T1, x2::T1, p0::T2)::Float64
    supp = _support(estimator.design)
    if !_isInSupport(supp, x1, x2)
        return(NaN)
    end
    estimates     = estimate.([estimator], supp[:, 1], supp[:, 2])
    estimateobs   = estimate(estimator, x1, x2)
    ismoreextreme = estimates .>= estimateobs
    return probability.([estimator.design], supp[ismoreextreme, 1], supp[ismoreextreme, 2], p0) |> sum
end

function bias{T<:Real}(estimator::BinaryTwoStageDesignEstimator, p::T)::Float64
    supp      = _support(estimator.design)
    estimates = estimate.([estimator], supp[:, 1], supp[:, 2])
    probs     = probability.([estimator.design], supp[:, 1], supp[:, 2], p)
    return dot(probs, estimates) - p
end

function rmse{T<:Real}(estimator::BinaryTwoStageDesignEstimator, p::T)::Float64
    supp   = _support(estimator.design)
    errors = estimate.([estimator], supp[:, 1], supp[:, 2]) .- p
    probs  = probability.([estimator.design], supp[:, 1], supp[:, 2], p)
    return sqrt(dot(probs, errors.^2))
end
