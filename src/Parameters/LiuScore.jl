type LiuScore{
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
    function LiuScore(
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
        (fs >= 0) & (fs < 1) ? nothing : throw(InexactError())
        (fp >= 0) & (fp < 1) ? nothing : throw(InexactError())
        new(samplespace, p0, pmcrv, prior, alpha, beta, mincondpower, fs, fp, smoothness)
    end
end
function LiuScore{T_samplespace<:SampleSpace}(
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
    LiuScore{T_samplespace, T_gs, T_eff, T_mcp}(
        samplespace, p0, pmcrv, prior, alpha, beta, minconditionalpower, fs, fp, smoothness
    )
end

maxsamplesize(params::LiuScore) = maxsamplesize(params.samplespace)

allowsstoppingforefficacy{T_samplespace, T_gs, T_eff, T_mcp}(params::LiuScore{T_samplespace, T_gs, T_eff, T_mcp}) = T_eff == AllowStoppingForEfficacy ? true : false

isgroupsequential{T_samplespace, T_gs, T_eff, T_mcp}(params::LiuScore{T_samplespace, T_gs, T_eff, T_mcp}) = T_gs == GroupSequential ? true : false

hasmonotoneconditionalpower{T_samplespace, T_gs, T_eff, T_mcp}(params::LiuScore{T_samplespace, T_gs, T_eff, T_mcp}) = T_mcp == MonotoneConditionalPower ? true : false

minconditionalpower(params::LiuScore) = params.mincondpower

beta(params::LiuScore) = params.beta

smoothness(params::LiuScore) = params.smoothness

function ros{T_P<:LiuScore}(design::AbstractBinaryTwoStageDesign, params::T_P, p1::Real)
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
    return min(100, ros/(1 - params.fs))
end
function rup{T_P<:LiuScore}(design::AbstractBinaryTwoStageDesign, params::T_P, p1::Real)
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

score{T_P<:LiuScore}(design::AbstractBinaryTwoStageDesign, params::T_P, p1::Real) = ros(design, params, p1) + rup(design, params, p1)

function score{T_P<:LiuScore}(design::AbstractBinaryTwoStageDesign, params::T_P)
    function f(p)
        res = params.prior(p)*score(design, params, p)
        return min(100, res)
    end
    return quadgk(f, params.pmcrv, 1, reltol = .001, maxevals = 1e5)[1]
end


function _createProblem{T<:Integer}(
    n1::T,      # stage one sample size
    params::LiuScore;
    npivots      = 15,
    npriorpivots = 25
)
    ss = samplespace(params)
    nmax = maxsamplesize(ss)
    possible(n1, ss) ? nothing : throw(InexactError())
    a  = alpha(params)
    b  = params.beta
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
    priorpivots = collect(linspace(params.pmcrv, 1.0, npriorpivots + 2))[2:(npriorpivots + 1)] # leave out boundary values!
    dp          = priorpivots[2] - priorpivots[1]
    priorvals   = prior.(priorpivots) ./ sum(prior.(priorpivots) .* dp) # normalize to 1
    # add type one error rate constraint
    @constraint(m,
        sum(dbinom(x1, n1, p0)*_cpr(x1, n1, n, c, p0)*y[x1, n, c] for
            x1 in 0:n1, n in nvals, c in cvals
        ) <= a
    )
    # add conditional type two error rate constraint (power)
    for x1 in 0:n1
        posterior1 = priorvals .* dbinom.(x1, n1, priorpivots) .* dp
        z1 = sum(posterior1)
        @constraint(m,
            sum(sum(posterior1 .* _cpr.(x1, n1, n, c, priorpivots))/z1*y[x1, n, c] for
                n in nvals, c in cvalsfinite
            ) + sum(y[x1, n, c] for
                n in nvals, c in cvalsinfinite
            ) >= minconditionalpower(params)
        )
        # ensure monotonicity if required
        if x1 >= 1 & hasmonotoneconditionalpower(params)
            posterior2 = priorvals .* dbinom.(x1 - 1, n1, priorpivots) .* dp
            z2 = sum(posterior2)
            @constraint(m,
                sum(sum(posterior1 .* _cpr.(x1, n1, n, c, priorpivots))/z1*y[x1, n, c]
                        - sum(posterior2 .* _cpr.(x1 - 1, n1, n, c, priorpivots))/z2*y[x1 - 1, n, c] for
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
            absdiff[x1] >= sum(n*(y[x1, n, c] - y[x1 + 1, n, c]) for n in nvals, c in cvalsfinite) - sum(2*nmax*y[x1, n, c] for n in nvals, c in cvalsinfinite) - sum(2*nmax*y[x1 + 1, n, c] for n in nvals, c in cvalsinfinite)
        )
        @constraint(m,
            -absdiff[x1] <= -sum(n*(y[x1, n, c] - y[x1 + 1, n, c]) for n in nvals, c in cvalsfinite) + sum(2*nmax*y[x1, n, c] for n in nvals, c in cvalsinfinite) + sum(2*nmax*y[x1 + 1, n, c] for n in nvals, c in cvalsinfinite)
        )
        @constraint(m, absdiff[x1] <= maxdiff)
    end
    # add optimality criterion
    nreq__(p, powerreq) = nreq_(p, powerreq, p0, a)
    # construct expressions for ROS, upper bound necessary to guarantee finteness!
    @expression(m, ros[p in priorpivots],
        1/(1 - params.fs) * sum(
            dbinom(x1, n1, p) * max(0, min(100.0, (n / nreq__(p, 1 - b) - 1))) * y[x1, n, c] for
            x1 in 0:n1, n in nvals, c in cvals
        )
    )
    # construct expressions for power
    @expression(m, designpower[p in priorpivots],
        sum(dbinom(x1, n1, p) * _cpr(x1, n1, n, c, p) * y[x1, n, c] for
            x1 in 0:n1, n in nvals, c in cvals
        )
    )
    pivots = collect(linspace(0, 1, npivots)) # lambda formulaion requires edges!
    @variable(m, 0 <= lambda[priorpivots, pivots] <= 1)
    function frup(power, p; rupmax = 100.0)
        nom   = max(0, nreq__(p, 1 - b) - nreq__(p, power))
        denom = nreq__(p, 1 - b) - nreq__(p, (1 - params.fp)*(1 - b))
        if (denom == 0.0) & (nom > 0)
            return rupmax
        end
        if (denom == 0.0) & (nom == 0)
            return 0.0
        end
        res = max(0, min(rupmax, nom / denom))
        return res
    end
    @variable(m, rup[priorpivots])
    for p in priorpivots
        addSOS2(m, [lambda[p, piv] for piv in pivots])
        @constraint(m, sum(lambda[p, piv] for piv in pivots) == 1)
        @constraint(m, sum(lambda[p, piv]*piv for piv in pivots) == designpower[p]) # defines lambdas!
        @constraint(m, sum(lambda[p, piv]*frup(piv, p) for piv in pivots) == rup[p])
    end
    @objective(m, Min,
        sum(priorvals[i]*(
            rup[priorpivots[i]] +
            ros[priorpivots[i]]
            )*dp for i in 1:npriorpivots) # normalized version
    )
    return m, y
end

function _isfeasible(design::BinaryTwoStageDesign, params::LiuScore)
    return true
end


# utility
function nreq_(p, power, p0, alpha; nmax = 10000.0)
    if alpha >= power
        return 0.0 # qnorm(1 - alpha) + qnorm(power) < 0 ....
    end
    if p <= p0
        return 0.0 # no effect
    end
    if power == 0
        return 0.0
    end
    if abs(power - 1.0) < 5*eps(Float64)
        return nmax
    end
    res = (1 - p)*p*( (qnorm(1 - alpha) + qnorm(power)) / (p0 - p) )^2
    res = min(nmax, max(5, res))
    return res
end
