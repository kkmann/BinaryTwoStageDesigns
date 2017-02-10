abstract BinaryTwoStageDesignEstimator

# these two enable array broadcasting of all methods!
Base.size(::BinaryTwoStageDesignEstimator) = ()
Base.getindex(estimator::BinaryTwoStageDesignEstimator, i) = estimator

design(estimator::BinaryTwoStageDesignEstimator) = try return estimator.design catch error("not implemented") end

estimate{T<:Integer}(estimator::BinaryTwoStageDesignEstimator, x1::T, x2::T) = error("not implemented")

function p{T1<:Integer, T2<:Real}(estimator::BinaryTwoStageDesignEstimator, x1::T1, x2::T1, p0::T2)
    checkx1x2(x1, x2, design(estimator))
    supp          = support(design(estimator))
    estimates     = estimate.(estimator, supp[:, 1], supp[:, 2])
    estimateobs   = estimate(estimator, x1, x2)
    ismoreextreme = estimates .>= estimateobs
    return pdf.(design(estimator), supp[ismoreextreme, 1], supp[ismoreextreme, 2], p0) |> sum
end

function bias{T<:Real}(estimator::BinaryTwoStageDesignEstimator, p::T)
    supp      = support(design(estimator))
    estimates = estimate.(estimator, supp[:, 1], supp[:, 2])
    probs     = pdf.(design(estimator), supp[:, 1], supp[:, 2], p)
    return dot(probs, estimates) - p
end

function rmse{T<:Real}(estimator::BinaryTwoStageDesignEstimator, p::T)
    supp     = support(design(estimator))
    return sqrt( sum(
        pdf.(design(estimator), supp[:, 1], supp[:, 2], p) .*
            (estimate.(estimator, supp[:, 1], supp[:, 2]) .- p).^2
    ) )
end

# ToDo: get inconsistencies
