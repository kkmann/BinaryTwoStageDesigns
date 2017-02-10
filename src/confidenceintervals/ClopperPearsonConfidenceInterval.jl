type ClopperPearsonConfidenceInterval <: ConfidenceInterval
    estimator::BinaryTwoStageDesignEstimator
    confidence::Float64

    function ClopperPearsonConfidenceInterval{T<:Real}(
        estimator::BinaryTwoStageDesignEstimator;
        confidence::T = .9
    )
        checkp(confidence)
        new(estimator, confidence)
    end
end

estimator(ci::ClopperPearsonConfidenceInterval) = ci.estimator

design(ci::ClopperPearsonConfidenceInterval) = ci |> estimator |> design


function limits{T<:Integer}(ci::ClopperPearsonConfidenceInterval, x1::T, x2::T; k::Integer = 1001)
    ispossible(design(ci), x1, x2) ? nothing : throw(InexactError())
    n1     = interimsamplesize(design(ci))
    nmax   = maximum(samplesize(design(ci)))
    grid   = collect(linspace(0, 1, k))
    supp   = support(design(ci)) # compute p values by hand is faster than calling p()
    tmp    = zeros(Float64, size(supp, 1), k)
    estimates = estimate.(estimator(ci), supp[:, 1], supp[:, 2])
    for i in 1:size(supp, 1)
        tmp[i, :] = pdf.(design(ci), supp[i, 1], supp[i, 2], linspace(0, 1, k))
    end
    probs = [supp samplesize.(design(ci), supp[:, 1]) estimates tmp]
    probs = sortrows(probs, by = x -> x[4], rev = true) # sort in descending order of estimates
    # compute p vals for H0: p < p0
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
    if x1 + x2 == samplesize(design(ci), x1)
        limits[2] = 1.0
    end
    return limits
end
