type EB{
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
    a
    b
    k
    smoothness
    function EB(
        samplespace,
        p0, pmcrv, prior,
        alpha,
        mincondpower,
        gamma, lambda,
        a, b,
        k,
        smoothness
    )
        any(!([alpha; p0; pmcrv; mincondpower] .>= 0.0)) ? throw(InexactError()) : nothing
        any(!([alpha; p0; pmcrv; mincondpower] .<= 1.0)) ? throw(InexactError()) : nothing
        abs(quadgk(prior, 0, 1)[1] - 1) <= .001 ? nothing: throw(InexactError())
        new(samplespace, p0, pmcrv, prior, alpha, mincondpower, gamma, lambda, a, b, k, smoothness)
    end
end
function EB{T_samplespace<:SampleSpace}(
    samplespace::T_samplespace,
    p0, pmcrv, prior, alpha,
    gamma, lambda;
    a::Real                        = 1,
    b::Real                        = 1,
    k::Real                        = 1,
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
    EB{T_samplespace, T_gs, T_eff, T_mcp}(
        samplespace, p0, pmcrv, prior, alpha, minconditionalpower, gamma, lambda, a, b, k, smoothness
    )
end

maxsamplesize(params::EB) = maxsamplesize(params.samplespace)

allowsstoppingforefficacy{T_samplespace, T_gs, T_eff, T_mcp}(params::EB{T_samplespace, T_gs, T_eff, T_mcp}) = T_eff == AllowStoppingForEfficacy ? true : false

isgroupsequential{T_samplespace, T_gs, T_eff, T_mcp}(params::EB{T_samplespace, T_gs, T_eff, T_mcp}) = T_gs == GroupSequential ? true : false

hasmonotoneconditionalpower{T_samplespace, T_gs, T_eff, T_mcp}(params::EB{T_samplespace, T_gs, T_eff, T_mcp}) = T_mcp == MonotoneConditionalPower ? true : false

minconditionalpower(params::EB) = params.mincondpower

smoothness(params::EB) = params.smoothness

function g(params::EB, power)
    power >= 0 ? nothing : error("power cannot be negative")
    power <= 1 ? nothing : error("power cannot be greater than 1")
    return Distributions.cdf(Distributions.Beta(params.a, params.b), power)
end

function score(design::AbstractBinaryTwoStageDesign, params::EB, p::Real)
    ess  = SampleSize(design, p) |> mean
    wpwr = p < mcrv(params) ? 0.0 : g(params, power(design, p))
    return -*(params.lambda*wpwr - params.gamma*(params.k*(1 - p) + p)*ess)
end

function score(design::AbstractBinaryTwoStageDesign, params::EB)
    return quadgk(p -> params.prior(p)*score(design, params, p), 0, 1, reltol = .0001)[1]
end

function _createProblem{T<:Integer}(
    n1::T,      # stage one sample size
    params::EB
)
    ss     = samplespace(params)
    nmax   = maxsamplesize(ss, n1)
    possible(n1, ss) ? nothing : throw(InexactError())
    a      = alpha(params)
    p0     = null(params)
    pmcrv  = params.pmcrv
    prior  = params.prior
    gamma  = params.gamma
    lambda = params.lambda
    maxdiff = min(2*nmax, smoothness(params)) # make finite!
    # define base problem
    m, y = _createBaseProblem(n1, params) # c.f. util.jl
    nvals = getnvals(ss, n1)
    cvalsfinite = 0:(maximum(nvals) - 1)
    if allowsstoppingforefficacy(params)
        cvals = [-Inf; cvalsfinite; Inf]
        cvalsinfinite = [-Inf; Inf]
    else
        cvals = [cvalsfinite; Inf]
        cvalsinfinite = [Inf]
    end
    @constraint(m, # type one error rate constraint
        sum(dbinom(x1, n1, p0)*_cpr(x1, n1, n, c, p0)*y[x1, n, c] for
            x1 in 0:n1, n in nvals, c in cvals
        ) <= alpha(params)
    )
    # precompute prior probabilities of X1=x1
    prior_prob = zeros(n1 + 1)
    response_prob = zeros(n1 + 1)
    for x1 in 0:n1
        prior_prob[x1 + 1]    = quadgk(p ->   prior(p)*dbinom(x1, n1, p), 0, 1, reltol = 0.0001)[1]
        response_prob[x1 + 1] = quadgk(p -> p*prior(p)*dbinom(x1, n1, p), 0, 1, reltol = 0.0001)[1] / prior_prob[x1 + 1]
    end
    @expression(m, ec, # expected cost
        sum(
            prior_prob[x1 + 1] * gamma * n * (params.k*(1 - response_prob[x1 + 1]) + response_prob[x1 + 1]) * y[x1, n , c]
            for x1 in 0:n1, n in nvals, c in cvals
        )
    )
    z = zeros(n1 + 1) # normalizing constants of the stage-one conditional power
    for x1 in 0:n1
        z[x1 + 1] = quadgk(p -> prior(p)*dbinom(x1, n1, p), pmcrv, 1, reltol = 0.0001)[1]
    end
    @expression(m, cewp[x1 in 0:n1], # conditional expected weighted power
        sum(
            quadgk(p ->
                    prior(p)*dbinom(x1, n1, p)/z[x1 + 1]* # stage one posterior f(p|x1, p >= pmcrv)
                    g(params, _cpr(x1, n1, n, c, p)), # weighted power
                pmcrv, 1, reltol = 0.0001
            )[1] * y[x1, n, c]
            for n in nvals, c in cvals
        )
    )
    for x1 in 0:n1
        @constraint(m, # add conditional type two error rate constraint (power)
            cewp[x1] >= minconditionalpower(params)
        )
        if x1 >= 1 & hasmonotoneconditionalpower(params)
            @constraint(m, # ensure monotonicity if required
                cewp[x1] - cewp[x1 - 1] >= 0
            )
        end
    end
    # pre compute prior probabilities of X1=x1 given p > pmcrv
    z = quadgk(p -> prior(p), pmcrv, 1, reltol = 0.0001)[1]
    prior_prob_pmcrv = zeros(n1 + 1)
    for x1 in 0:n1
        prior_prob_pmcrv[x1 + 1] = quadgk(p -> prior(p)*dbinom(x1, n1, p), pmcrv, 1, reltol = 0.0001)[1] / z
    end
    @expression(m, ewp, # expected weighted power
        sum(
            prior_prob_pmcrv[x1 + 1] * cewp[x1]
            for x1 in 0:n1
        )
    )
    @objective(m, Min,
        -(params.lambda * quadgk(p -> prior(p), pmcrv, 1, reltol = 0.0001)[1] * ewp - ec)
    )
    return m, y
end

function _isfeasible(design::BinaryTwoStageDesign, params::EB)
    return true
end
