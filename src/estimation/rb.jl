type RaoBlackwellizedEstimator <: BinaryTwoStageDesignEstimator
    design::AbstractBinaryTwoStageDesign

    function RaoBlackwellizedEstimator(
        design::AbstractBinaryTwoStageDesign
    )
        new(design)
    end
end

function estimate{T<:Integer}(estimator::RaoBlackwellizedEstimator, x1::T, x2::T)::Float64
    supp = _support(estimator.design)
    if !_isInSupport(supp, x1, x2)
        return NaN
    else
        n1   = getInterimSampleSize(estimator.design)
        nmax = maximum(getSampleSize(estimator.design))
        nominator   = BigInt(0.0) # we need big ints, as binomials get HUGE!
        denominator = BigInt(0.0)
        x = x1 + x2
        n = getSampleSize(estimator.design, x1)
        for xx1 in 0:(min(x, n1))
            if getSampleSize(estimator.design, xx1) == n
                nominator   += BigInt(binomial(n1 - 1, xx1 - 1))
                denominator += BigInt(binomial(n1, xx1))
            end
        end
        return nominator/denominator
    end
end
