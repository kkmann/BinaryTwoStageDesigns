type BMESS{
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
    beta
    a
    b
    mincondpower
    smoothness
    function BMESS(
        samplespace,
        p0, pmcrv, prior,
        alpha, beta,
        a, b,
        mincondpower, smoothness
    )
        alpha >= 0 ? nothing : error("alpha must be >= 0")
        alpha <= 1 ? nothing : error("alpha must be <= 1")
        beta  >= 0 ? nothing : error("beta must be >= 0")
        beta  <= 1 ? nothing : error("beta must be <= 1")
        p0    >= 0 ? nothing : error("p0 must be >= 0")
        p0    <= 1 ? nothing : error("p0 must be <= 1")
        pmcrv >= p0 ? nothing : error("pmcrv must be >= p0")
        pmcrv <= 1 ? nothing : error("pmcrv must be <= 1")
        a     >= 1 ? nothing : error("a must be >= 1")
        b     >= 1 ? nothing : error("b must be >= 1")
        mincondpower >= 0 ? nothing : error("mincondpower must be >= 0")
        mincondpower <= 1 ? nothing : error("mincondpower must be <= 1")
        abs(quadgk(prior, 0, 1)[1] - 1) <= .001 ? nothing: error("prior should integrate to 1")
        new(samplespace, p0, pmcrv, prior, alpha, beta, a, b, mincondpower, smoothness)
    end
end
function BMESS{T_samplespace<:SampleSpace}(
    samplespace::T_samplespace,
    p0, pmcrv, prior,
    alpha, beta;
    a::Real                        = 1,
    b::Real                        = 1,
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
    BMESS{T_samplespace, T_gs, T_eff, T_mcp}(
        samplespace, p0, pmcrv, prior, alpha, beta, a, b, minconditionalpower, smoothness
    )
end

maxsamplesize(params::BMESS) = maxsamplesize(params.samplespace)

allowsstoppingforefficacy{T_samplespace, T_gs, T_eff, T_mcp}(params::BMESS{T_samplespace, T_gs, T_eff, T_mcp}) = T_eff == AllowStoppingForEfficacy ? true : false

isgroupsequential{T_samplespace, T_gs, T_eff, T_mcp}(params::BMESS{T_samplespace, T_gs, T_eff, T_mcp}) = T_gs == GroupSequential ? true : false

hasmonotoneconditionalpower{T_samplespace, T_gs, T_eff, T_mcp}(params::BMESS{T_samplespace, T_gs, T_eff, T_mcp}) = T_mcp == MonotoneConditionalPower ? true : false

minconditionalpower(params::BMESS) = params.mincondpower

smoothness(params::BMESS) = params.smoothness

function score(design::AbstractBinaryTwoStageDesign, params::BMESS, p::Real)
    ss = SampleSize(design, p)
    return mean(ss)
end
function score(design::AbstractBinaryTwoStageDesign, params::BMESS)
    return quadgk(p -> prior(params, p)*score(design, params, p), 0, 1, reltol = .001)[1]
end

function _createProblem{T<:Integer}(
    n1::T,      # stage one sample size
    params::BMESS
)
    ss = samplespace(params)
    possible(n1, ss) ? nothing : error("n1 not compatible with sample space")
    nmax    = maxsamplesize(ss, n1)
    p0      = null(params)
    pmcrv   = params.pmcrv
    prior   = params.prior
    maxdiff = min(2*nmax, smoothness(params)) # make finite!
    # define base problem
    m, y    = _createBaseProblem(n1, params) # c.f. util.jl # TODO: implement smoothness here!
    nvals   = getnvals(ss, n1)
    cvalsfinite = 0:(maximum(nvals) - 1)
    if allowsstoppingforefficacy(params)
        cvals = [-Inf; cvalsfinite; Inf]
        cvalsinfinite = [-Inf; Inf]
    else
        cvals = [cvalsfinite; Inf]
        cvalsinfinite = [Inf]
    end
    g(pow) = Distributions.cdf(Distributions.Beta(params.a, params.b), pow)

    @constraint(m, # add type one error rate constraint
        sum(dbinom(x1, n1, p0)*_cpr(x1, n1, n, c, p0)*y[x1, n, c] for
            x1 in 0:n1, n in nvals, c in cvals
        ) <= alpha(params)
    )
    @expression(m, cewp[x1 in 0:n1], # define expressions for conditional weighted power given X1=x1
        sum(quadgk(p -> prior(p)*dbinom(x1, n1, p)*g(_cpr.(x1, n1, n, c, p)), pmcrv, 1, reltol = 0.001)[1] /
                quadgk(p -> prior(p)*dbinom(x1, n1, p), pmcrv, 1, reltol = 0.001)[1] * y[x1, n, c] for
            n in nvals, c in cvals
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
    @constraint(m,
        sum(cewp[x1] * quadgk(p -> dbinom(x1, n1, p)*prior(p), 0, 1, reltol = 0.001)[1] for
            x1 in 0:n1
        ) >= 1 - params.beta
    )
    # # add constraints for smoothness (maxmial difference between successive
    # # values of n on the continuation set)
    # @variable(m, absdiff[x1 = 0:(n1 - 1)])
    # for x1 in 0:(n1 - 1)
    #     if allowsstoppingforefficacy(params)
    #         tmp = [-Inf; cvalsfinite]
    #     else
    #         tmp = cvalsfinite
    #     end
    #     @constraint(m,
    #         absdiff[x1] >= sum(n*(y[x1, n, c] - y[x1 + 1, n, c]) for n in nvals, c in tmp)
    #             - sum(2*nmax*y[x1 + i, n, Inf] for n in nvals, i in 0:1) # disables constraint for stopping for futility
    #     )
    #     @constraint(m,
    #         -absdiff[x1] <= -sum(n*(y[x1, n, c] - y[x1 + 1, n, c]) for n in nvals, c in tmp)
    #             + sum(2*nmax*y[x1 + 1, n, Inf] for n in nvals, i in 0:1)
    #     )
    #     @constraint(m, absdiff[x1] <= maxdiff)
    # end
    # add optimality criterion
    # expected costs
    @objective(m, Min,
        sum(n * quadgk(p -> dbinom(x1, n1, p)*prior(p), 0, 1, reltol = 0.001)[1] * y[x1, n, c] for # (1 - Distributions.cdf(Distributions.Beta(1, 6.475), p))*
            x1 in 0:n1, n in nvals, c in cvals
        )
    )
    return m, y
end

function _isfeasible(design::BinaryTwoStageDesign, params::BMESS)
    return true
end
