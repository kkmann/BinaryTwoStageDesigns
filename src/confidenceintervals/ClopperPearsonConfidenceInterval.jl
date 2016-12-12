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
    # it is faster to compute the pvalues by hand...
    supp = BinaryTwoStageDesigns._support(getDesign(ci))
    tmp = zeros(Float64, size(supp, 1), k)
    estimates = zeros(Float64, size(supp, 1))
    prange = linspace(0, 1, k)
    for i in 1:size(supp, 1)
        estimates[i] = estimate(ci.estimator, supp[i, 1], supp[i, 2])
        for j in 1:k
            tmp[i, j] = probability(getDesign(ci), supp[i, 1], supp[i, 2], prange[j])
        end
    end
    probs = [supp getSampleSize.(getDesign(ci), supp[:, 1]) estimates tmp]
    probs = sortrows(probs, by = x -> x[4], rev = true) # sort in descending order of estimates
    # compute pvals for H0: p < p0
    pvals_leq = probs |> x -> cumsum(x, 1)
    pvals_leq[:, 1:4] = probs[:, 1:4]
    pvals_leq = pvals_leq[(pvals_leq[:, 1] .== x1) & (pvals_leq[:, 2] .== x2), 5:(k + 4)][1, :]

    # compute pvals for H0: p > p0
    pvals_geq = sortrows(probs, by = x -> x[4]) |> x -> cumsum(x, 1)
    pvals_geq[:, 1:4] = sortrows(probs, by = x -> x[4])[:, 1:4]
    pvals_geq = sortrows(pvals_geq, by = x -> x[4], rev = true)
    pvals_geq = pvals_geq[(pvals_geq[:, 1] .== x1) & (pvals_geq[:, 2] .== x2), 5:(k + 4)][1, :]

    indlow = findfirst(pvals_leq .>= (1 - ci.confidence)/2.0)
    indlow = indlow == 0 ? 1 : indlow
    indup  = findlast(pvals_geq .>= (1 - ci.confidence)/2.0)
    indup  = indup == 0 ? k : indup
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

function coverage{T<:Integer}(
    ci::BinaryTwoStageDesignConfidenceInterval; k::T = 1001
)
    design = getDesign(ci)
    n1     = getInterimSampleSize(design)
    nmax   = maximum(getSampleSize(design))
    # construct array of possible outcomes
    supp   = _support(design)
    supp = BinaryTwoStageDesigns._support(getDesign(ci))
    tmp = zeros(Float64, size(supp, 1), k)
    estimates = zeros(Float64, size(supp, 1))
    prange = linspace(0, 1, k)
    for i in 1:size(supp, 1)
        estimates[i] = estimate(ci.estimator, supp[i, 1], supp[i, 2])
        for j in 1:k
            tmp[i, j] = probability(getDesign(ci), supp[i, 1], supp[i, 2], prange[j])
        end
    end
    probs = [supp getSampleSize.(getDesign(ci), supp[:, 1]) estimates tmp]
    probs = sortrows(probs, by = x -> x[4], rev = true) # sort in descending order of estimates
    # compute pvals for H0: p < p0
    pvals_leq = probs |> x -> cumsum(x, 1)
    pvals_leq[:, 1:4] = probs[:, 1:4]
    # compute pvals for H0: p > p0
    pvals_geq = sortrows(probs, by = x -> x[4]) |> x -> cumsum(x, 1)
    pvals_geq[:, 1:4] = sortrows(probs, by = x -> x[4])[:, 1:4]
    pvals_geq = sortrows(pvals_geq, by = x -> x[4], rev = true)
    #
    tmp2 = (
        ((pvals_leq[:, 5:(k + 4)] .>= .05) & (pvals_geq[:, 5:(k + 4)] .>= .05)) |
        ((pvals_leq[:, 5:(k + 4)] .>= .05) & repmat(sum(pvals_leq[:, 1:2], 2) .== 0, 1, k)) |
        ((pvals_geq[:, 5:(k + 4)] .>= .05) & repmat(sum(pvals_geq[:, 1:2], 2) .== pvals_geq[:, 3], 1, k))
    )  .* probs[:, 5:(k + 4)] |> x -> sum(x, 1)
    return DataFrame(presp = linspace(0, 1, k), coverage = tmp2[1, :])
end

function getEstimator(ci::ClopperPearsonConfidenceInterval)
    return ci.estimator
end

function getDesign(ci::ClopperPearsonConfidenceInterval)
    return ci.estimator.design
end
