type Liu{
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
    mincondpower
    fs
    fp
    smoothness
    function Liu(
        samplespace,
        p0, pmcrv, prior,
        alpha, beta,
        mincondpower,
        fs, fp,
        smoothness
    )
        any(!([alpha; beta; p0; pmcrv; mincondpower] .>= 0.0)) ? throw(InexactError()) : nothing
        any(!([alpha; beta; p0; pmcrv; mincondpower] .<= 1.0)) ? throw(InexactError()) : nothing
        quadgk(prior, 0, pmcrv)[1] <= 0.001 ? nothing: throw(InexactError())
        quadgk(prior, pmcrv, 1)[1] >= 0.999 ? nothing: throw(InexactError())
        (fp >= 0) & (fp < 1) ? nothing : throw(InexactError())
        fs > 1 ? nothing : throw(InexactError())
        new(samplespace, p0, pmcrv, prior, alpha, beta, mincondpower, fs, fp, smoothness)
    end
end
function Liu{T_samplespace<:SampleSpace}(
    samplespace::T_samplespace,
    p0, pmcrv, prior, alpha, beta,
    fs, fp;
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
    Liu{T_samplespace, T_gs, T_eff, T_mcp}(
        samplespace, p0, pmcrv, prior, alpha, beta, minconditionalpower, fs, fp, smoothness
    )
end

maxsamplesize(params::Liu) = maxsamplesize(params.samplespace)

allowsstoppingforefficacy{T_samplespace, T_gs, T_eff, T_mcp}(params::Liu{T_samplespace, T_gs, T_eff, T_mcp}) = T_eff == AllowStoppingForEfficacy ? true : false

isgroupsequential{T_samplespace, T_gs, T_eff, T_mcp}(params::Liu{T_samplespace, T_gs, T_eff, T_mcp}) = T_gs == GroupSequential ? true : false

hasmonotoneconditionalpower{T_samplespace, T_gs, T_eff, T_mcp}(params::Liu{T_samplespace, T_gs, T_eff, T_mcp}) = T_mcp == MonotoneConditionalPower ? true : false

minconditionalpower(params::Liu) = params.mincondpower

beta(params::Liu) = params.beta

smoothness(params::Liu) = params.smoothness

function ros{T_P<:Liu}(design::AbstractBinaryTwoStageDesign, params::T_P, p1::Real)
    powerreq = 1 - beta(params)
    n1 = interimsamplesize(design)
    nreq__(p, powerreq) = nreq_(p, powerreq, params.p0, params.alpha)
    denom = nreq__(p1, powerreq)
    if denom == 0
        return 100
    end
    ros = 0.0
    for x1 in 0:n1
        ros += dbinom(x1, n1, p1)*max(0.0, samplesize(design, x1)/denom - 1)
    end
    return min(100, ros/(params.fs - 1))
end
function rup{T_P<:Liu}(design::AbstractBinaryTwoStageDesign, params::T_P, p1::Real)
    if p1 < params.pmcrv
        return 0.0
    end
    powerreq = 1 - beta(params) # we can minimize this if we set power to zero!
    n1 = interimsamplesize(design)
    nreq__(p, powerreq) = nreq_(p, powerreq, params.p0, params.alpha)
    nom   = max(0, nreq__(p1, powerreq) - nreq__(p1, power(design, p1)))
    denom = nreq__(p1, powerreq) - nreq__(p1, (1 - params.fp)*powerreq)
    if (denom == 0.0) & (nom > 0)
        return 100
    end
    if (denom == 0.0) & (nom == 0)
        return 0.0
    end
    res = nom / denom
    return min(100, res)
end

score{T_P<:Liu}(design::AbstractBinaryTwoStageDesign, params::T_P, p1::Real) = ros(design, params, p1) + rup(design, params, p1)

function score{T_P<:Liu}(design::AbstractBinaryTwoStageDesign, params::T_P)
    function f(p)
        res = params.prior(p)*score(design, params, p)
        return min(100, res)
    end
    return quadgk(f, params.pmcrv, 1, reltol = .001)[1]
end


function _createProblem{T<:Integer}(
    n1::T,      # stage one sample size
    params::Liu
)
    ss      = samplespace(params)
    nmax    = maxsamplesize(ss, n1)
    possible(n1, ss) ? nothing : throw(InexactError())
    p0      = null(params)
    pmcrv   = mcrv(params)
    prior   = params.prior
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
    @expression(m, cep[x1 in 0:n1], # define expressions for conditional power given X1=x1
        sum(quadgk(p -> prior(p)*dbinom(x1, n1, p)*_cpr.(x1, n1, n, c, p), pmcrv, 1, reltol = 0.001)[1] /
                quadgk(p -> prior(p)*dbinom(x1, n1, p), pmcrv, 1, reltol = 0.001)[1] * y[x1, n, c] for
            n in nvals, c in cvals
        )
    )
    for x1 in 0:n1
        @constraint(m, # add conditional type two error rate constraint (power)
            cep[x1] >= minconditionalpower(params)
        )
        if x1 >= 1 & hasmonotoneconditionalpower(params)
            @constraint(m, # ensure monotonicity if required
                cep[x1] - cep[x1 - 1] >= 0
            )
        end
    end
    # precompute prior probabilities of X1=x1 given p >= pmcrv
    z = zeros(n1 + 1) # redefine as normalizing constants of the stage-one posterior
    prior_prob_pmcrv = zeros(n1 + 1)
    for x1 in 0:n1
        z[x1 + 1] = quadgk(p -> prior(p)*dbinom(x1, n1, p), pmcrv, 1, reltol = 0.001)[1]
        prior_prob_pmcrv[x1 + 1] = quadgk(p -> prior(p)*dbinom(x1, n1, p), pmcrv, 1, reltol = 0.001)[1]
    end
    nreq__(p, powerreq) = nreq_(p, powerreq, p0, alpha(params))
    function frup(power, p; rupmax = 100.0)
        if p < params.pmcrv
            return 0.0
        end
        nom   = max(0, nreq__(p, 1 - beta(params)) - nreq__(p, power))
        denom = nreq__(p, 1 - beta(params)) - nreq__(p, (1 - params.fp)*(1 - beta(params)))
        if (denom == 0.0) & (nom > 0)
            return rupmax
        end
        if (denom == 0.0) & (nom == 0)
            return 0.0
        end
        res = max(0, min(rupmax, nom / denom))
        return res
    end
    @expression(m, cros[x1 in 0:n1],
        1/(params.fs - 1) * sum(
            quadgk(p ->
                    prior(p)*dbinom(x1, n1, p)/z[x1 + 1]*
                    max(0, min(100.0, (n / nreq__(p, 1 - beta(params)) - 1))),
                pmcrv, 1, reltol = 0.001
            )[1]*y[x1, n, c]
            for n in nvals, c in cvals
        )
    )
    @expression(m, crup[x1 in 0:n1],
        sum(
            quadgk(p ->
                    prior(p)*dbinom(x1, n1, p)/z[x1 + 1]*
                    frup(_cpr(x1, n1, n, c, p), p),
                pmcrv, 1, reltol = 0.001
            )[1]*y[x1, n, c]
            for n in nvals, c in cvals
        )
    )
    @expression(m, ros, # expected ros
        sum(
            prior_prob_pmcrv[x1 + 1] * cros[x1]
            for x1 in 0:n1
        )
    )
    @expression(m, rup, # expected ros
        sum(
            prior_prob_pmcrv[x1 + 1] * crup[x1]
            for x1 in 0:n1
        )
    )
    @objective(m, Min,
        rup + ros
    )
    return m, y
end

function _isfeasible(design::BinaryTwoStageDesign, params::Liu)
    return true
end


# # utility
# function nreq_(p, power, p0, alpha; nmax = 10000.0)
#     if alpha >= power
#         return 0.0 # qnorm(1 - alpha) + qnorm(power) < 0 ....
#     end
#     if p <= p0
#         return 0.0 # no effect
#     end
#     if power == 0
#         return 0.0
#     end
#     if abs(power - 1.0) < 5*eps(Float64)
#         return nmax
#     end
#     res = (1 - p)*p*( (qnorm(1 - alpha) + qnorm(power)) / (p0 - p) )^2
#     res = min(nmax, max(0, res))
#     if isnan(res)
#         println(p, power, p0, alpha, nmax)
#     end
#     return res
# end
