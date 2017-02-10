type RaoBlackwellizedEstimator <: BinaryTwoStageDesignEstimator
    design::AbstractBinaryTwoStageDesign

    function RaoBlackwellizedEstimator(
        design::AbstractBinaryTwoStageDesign
    )
        new(design)
    end
end

function estimate{T<:Integer}(estimator::RaoBlackwellizedEstimator, x1::T, x2::T)
    checkx1x2(x1, x2, design(estimator))
    n1          = interimsamplesize(design(estimator))
    nominator   = BigInt(0.0) # we need big ints, as binomials get HUGE!
    denominator = BigInt(0.0)
    x = x1 + x2
    n = samplesize(design(estimator), x1)
    for xx1 in 0:(min(x, n1))
        if samplesize(design(estimator), xx1) == n
            nominator   += BigInt(binomial(n1 - 1, xx1 - 1))*BigInt(binomial(n - n1, x - xx1))
            denominator += BigInt(binomial(n1, xx1))*BigInt(binomial(n - n1, x - xx1))
        end
    end
    return convert(Float64, nominator/denominator)
end
