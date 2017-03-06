type KunzmannScore{
        T_samplespace<:SampleSpace,
        T_gs<:IsGroupSequential,
        T_eff<:StoppingForEfficacy,
        T_mcp<:HasMonotoneConditionalPower
} <: VagueAlternative
    samplespace::T_samplespace
    p0
    pmcrv
    prior
    alpha
    mincondpower
    gamma
    lambda
    minpow
    maxpow
    smoothness
    riskpremium
    function KunzmannScore(
        samplespace,
        p0, pmcrv, prior,
        alpha,
        mincondpower,
        gamma, lambda,
        minpow, maxpow,
        smoothness,
        riskpremium
    )
        any(!([alpha; minpow; maxpow; p0; pmcrv; mincondpower] .>= 0.0)) ? throw(InexactError()) : nothing
        any(!([alpha; minpow; maxpow; p0; pmcrv; mincondpower] .<= 1.0)) ? throw(InexactError()) : nothing
        abs(quadgk(prior, 0, 1)[1] - 1) <= .001 ? nothing: throw(InexactError())
        minpow >= maxpow ? throw(InexactError()) : nothing
        new(samplespace, p0, pmcrv, prior, alpha, mincondpower, gamma, lambda, minpow, maxpow, smoothness, riskpremium)
    end
end
function KunzmannScore{T_samplespace<:SampleSpace}(
    samplespace::T_samplespace,
    p0, pmcrv, prior, alpha,
    gamma, lambda;
    riskpremium::Real              = 0,
    minpow::Real                   = .7,
    maxpow::Real                   = .9,
    smoothness::Real               = Inf,
    minconditionalpower::Real      = 0.0,
    GROUPSEQUENTIAL::Bool          = false,
    STOPPINGFOREFFICACY::Bool      = true,
    MONOTONECONDITIONALPOWER::Bool = true
)
    T_gs = NotGroupSequential
    if GROUPSEQUENTIAL
        T_gs = GroupSequential
    end
    T_eff = AllowStoppingForEfficacy
    if !STOPPINGFOREFFICACY
        T_eff = NoStoppingForEfficacy
    end
    T_mcp = NoMonotoneConditionalPower
    if MONOTONECONDITIONALPOWER
        T_mcp = MonotoneConditionalPower
    end
    KunzmannScore{T_samplespace, T_gs, T_eff, T_mcp}(
        samplespace, p0, pmcrv, prior, alpha, minconditionalpower, gamma, lambda, minpow, maxpow, smoothness, riskpremium
    )
end

maxsamplesize(params::KunzmannScore) = maxsamplesize(params.samplespace)

allowsstoppingforefficacy{T_samplespace, T_gs, T_eff, T_mcp}(params::KunzmannScore{T_samplespace, T_gs, T_eff, T_mcp}) = T_eff == AllowStoppingForEfficacy ? true : false

isgroupsequential{T_samplespace, T_gs, T_eff, T_mcp}(params::KunzmannScore{T_samplespace, T_gs, T_eff, T_mcp}) = T_gs == GroupSequential ? true : false

hasmonotoneconditionalpower{T_samplespace, T_gs, T_eff, T_mcp}(params::KunzmannScore{T_samplespace, T_gs, T_eff, T_mcp}) = T_mcp == MonotoneConditionalPower ? true : false

minconditionalpower(params::KunzmannScore) = params.mincondpower

smoothness(params::KunzmannScore) = params.smoothness

function expectedcost(design::AbstractBinaryTwoStageDesign, params::KunzmannScore, p::Real)
    ss = SampleSize(design, p)
    return params.gamma * p * mean(ss) + (params.gamma + params.riskpremium) * (1 - p) * mean(ss)
end
function underpowerpenalty(design::AbstractBinaryTwoStageDesign, params::KunzmannScore, p::Real)
    if p < params.pmcrv
        return 0.0
    else
        return params.lambda * (1 - g(params, power(design, p)))
    end
end

score(design::AbstractBinaryTwoStageDesign, params::KunzmannScore, p::Real) = expectedcost(design, params, p) + underpowerpenalty(design, params, p)

function score(design::AbstractBinaryTwoStageDesign, params::KunzmannScore)
    f(p) = params.prior(p)*score(design, params, p)
    return quadgk(f, 0, 1, reltol = .001, maxevals = 1e5)[1]
end


function _createProblem{T<:Integer}(
    n1::T,      # stage one sample size
    params::KunzmannScore;
    npivots      = 15,
    npriorpivots = 50
)
    ss = samplespace(params)
    nmax = maxsamplesize(ss)
    possible(n1, ss) ? nothing : throw(InexactError())
    a  = alpha(params)
    p0 = null(params)
    prior = params.prior
    maxdiff = min(2*nmax, smoothness(params)) # make finite!
    # define base problem
    m, y = _createBaseProblem(n1, params) # c.f. util.jl
    nvals = n1:nmax
    cvalsfinite = 0:(nmax - 1)
    if allowsstoppingforefficacy(params)
        cvals         = [-Inf; cvalsfinite; Inf]
        cvalsinfinite = [-Inf; Inf]
    else
        cvals         = [cvalsfinite; Inf]
        cvalsinfinite = [Inf]
    end
    priorpivots  = collect(linspace(0.0, 1.0, npriorpivots + 2))[2:(npriorpivots + 1)] # leave out boundary values!
    dp           = priorpivots[2] - priorpivots[1]
    priorvals    = prior.(priorpivots) ./ sum(prior.(priorpivots) .* dp) # normalize to 1
    cpriorpivots = priorpivots[priorpivots .>= params.pmcrv]
    cpriorvals   = prior.(cpriorpivots) ./ sum(prior.(cpriorpivots) .* dp) # normalize to 1, dp stays the same
    # add type one error rate constraint
    @constraint(m,
        sum(dbinom(x1, n1, p0)*_cpr(x1, n1, n, c, p0)*y[x1, n, c] for
            x1 in 0:n1, n in nvals, c in cvals
        ) <= a
    )
    # add conditional type two error rate constraint (power)
    for x1 in 0:n1
        posterior1 = cpriorvals .* dbinom.(x1, n1, cpriorpivots) .* dp # this must be conditional on $\rho>\rho_mcrv
        z1 = sum(posterior1)
        @constraint(m,
            sum(sum(posterior1 .* _cpr.(x1, n1, n, c, cpriorpivots))/z1*y[x1, n, c] for
                n in nvals, c in cvalsfinite
            ) + sum(y[x1, n, c] for
                n in nvals, c in cvalsinfinite
            ) >= minconditionalpower(params)
        )
        # ensure monotonicity if required
        if x1 >= 1 & hasmonotoneconditionalpower(params)
            posterior2 = cpriorvals .* dbinom.(x1 - 1, n1, cpriorpivots) .* dp
            z2 = sum(posterior2)
            @constraint(m,
                sum(sum(posterior1 .* _cpr.(x1, n1, n, c, cpriorpivots))/z1*y[x1, n, c]
                        - sum(posterior2 .* _cpr.(x1 - 1, n1, n, c, cpriorpivots))/z2*y[x1 - 1, n, c] for
                    n in nvals, c in cvals
                ) >= 0
            )
        end
    end
    # add constraints for smoothness (maxmial difference between successive
    # values of n on the continuation set)
    @variable(m, absdiff[x1 = 0:(n1 - 1)])
    for x1 in 0:(n1 - 1)
        @constraint(m,
            absdiff[x1] >= sum(n*(y[x1, n, c] - y[x1 + 1, n, c]) for n in nvals, c in cvalsfinite)
            - sum(2*nmax*y[x1 + i, n, Inf] for n in nvals, i in 0:1) # disables constraint for stopping for futility
        )
        @constraint(m,
            -absdiff[x1] <= -sum(n*(y[x1, n, c] - y[x1 + 1, n, c]) for n in nvals, c in cvalsfinite) + sum(2*nmax*y[x1 + 1, n, Inf] for n in nvals, i in 0:1)
        )
        @constraint(m, absdiff[x1] <= maxdiff)
    end
    # add optimality criterion
    # expected costs
    @expression(m, ec[p in priorpivots],
        sum(
            (params.gamma * p + (params.gamma + params.riskpremium) * (1 - p)) * n * dbinom(x1, n1, p) * y[x1, n, c] for
            x1 in 0:n1, n in nvals, c in cvals
        )
    )
    # construct expressions for power
    @expression(m, designpower[p in priorpivots],
        sum(dbinom(x1, n1, p) * _cpr(x1, n1, n, c, p) * y[x1, n, c] for
            x1 in 0:n1, n in nvals, c in cvals
        )
    )
    pivots = [0; collect(linspace(params.minpow, params.maxpow, max(1, npivots - 2))); 1] # lambda formulaion requires edges! exploitn fact that function is locally linear!
    @variable(m, 0 <= lambdaSOS2[priorpivots, pivots] <= 1)
    @variable(m, upp[priorpivots]) # underpower penalty
    for p in priorpivots
        addSOS2(m, [lambdaSOS2[p, piv] for piv in pivots])
        @constraint(m, sum(lambdaSOS2[p, piv] for piv in pivots) == 1)
        @constraint(m, sum(lambdaSOS2[p, piv]*piv for piv in pivots) == designpower[p]) # defines lambdas!
        if p >= params.pmcrv
            @constraint(m, params.lambda * (1 - sum(lambdaSOS2[p, piv]*g(params, piv) for piv in pivots)) == upp[p])
        else
            @constraint(m, upp[p] == 0)
        end
    end
    @objective(m, Min,
        sum(priorvals[i]*(
            ec[priorpivots[i]] +
            upp[priorpivots[i]] #; priorpivots[i] >= params.pmcrv] # todo: check whether it is really that simple
            )*dp for i in 1:npriorpivots) # normalized version
    )
    return m, y
end

function _isfeasible(design::BinaryTwoStageDesign, params::KunzmannScore)
    return true
end

# utility
function g(params::KunzmannScore, power)
    power
    power < 0 ? throw(InexactError()) : nothing
    power > 1 ? throw(InexactError()) : nothing
    if power < params.minpow
        return power
    end
    if (power >= params.minpow) & (power <= params.maxpow)
        return params.minpow + (.5 + (params.maxpow - power)/(params.maxpow - params.minpow))*(power - params.minpow)
    end
    if (power > params.maxpow)
        return 1.0
    end
end
