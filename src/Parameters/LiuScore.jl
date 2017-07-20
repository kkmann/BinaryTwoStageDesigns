type LiuScore{T_samplespace<:SampleSpace} <: VagueAlternative
    samplespace::T_samplespace
    p0
    pmcrv
    prior
    alpha
    beta
    mincondpower
    fs
    fp
    MONOTONECONDITIONALPOWER
    npriorpivots # number of pivots for prior evaluation
    npivots # number of pivots for approximation of g
    function LiuScore(
        samplespace,
        p0, pmcrv, prior,
        alpha, beta,
        mincondpower,
        fs, fp,
        MONOTONECONDITIONALPOWER,
        npriorpivots, # number of pivots for prior evaluation
        npivots, # number of pivots for approximation of g
    )
        any(!([alpha; beta; p0; pmcrv; mincondpower] .>= 0.0)) ? throw(InexactError()) : nothing
        any(!([alpha; beta; p0; pmcrv; mincondpower] .<= 1.0)) ? throw(InexactError()) : nothing
        quadgk(prior, 0, pmcrv)[1] <= 0.001 ? nothing: throw(InexactError())
        quadgk(prior, pmcrv, 1)[1] >= 0.999 ? nothing: throw(InexactError())
        (fp >= 0) & (fp < 1) ? nothing : throw(InexactError())
        fs > 1 ? nothing : throw(InexactError())
        new(samplespace, p0, pmcrv, prior, alpha, beta, mincondpower, fs, fp,
            MONOTONECONDITIONALPOWER, npriorpivots, npivots
        )
    end
end
function LiuScore{T_samplespace<:SampleSpace}(
    samplespace::T_samplespace,
    p0, pmcrv, prior, alpha, beta,
    fs, fp;
    npriorpivots = 50,
    npivots = 15,
    minconditionalpower::Real = 0.0,
    MONOTONECONDITIONALPOWER::Bool = false,
)
    LiuScore{T_samplespace}(
        samplespace, p0, pmcrv, prior, alpha, beta, minconditionalpower, fs, fp, MONOTONECONDITIONALPOWER, npriorpivots, npivots
    )
end

maxsamplesize(params::LiuScore) = maxsamplesize(params.samplespace)

hasmonotoneconditionalpower(params::LiuScore) = params.MONOTONECONDITIONALPOWER

minconditionalpower(params::LiuScore) = params.mincondpower

beta(params::LiuScore) = params.beta

function ros{T_P<:LiuScore}(design::BinaryTwoStageDesign, params::T_P, p1::Real)
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
function rup{T_P<:LiuScore}(design::BinaryTwoStageDesign, params::T_P, p1::Real)
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

score{T_P<:LiuScore}(design::BinaryTwoStageDesign, params::T_P, p1::Real) = ros(design, params, p1) + rup(design, params, p1)

function score{T_P<:LiuScore}(design::BinaryTwoStageDesign, params::T_P)
    function f(p)
        res = params.prior(p)*score(design, params, p)
        return min(100, res)
    end
    return quadgk(f, params.pmcrv, 1, reltol = .001, maxevals = 1e5)[1]
end

function completemodel{T<:Integer}(ipm::IPModel, params::LiuScore, n1::T)
    ss = samplespace(params)
    !possible(n1, ss) ? error("n1 and sample space incompatible") : nothing
    # extract ip model
    m             = ipm.m
    y             = ipm.y
    nvals         = ipm.nvals
    cvals         = ipm.cvals
    cvalsfinite   = ipm.cvalsfinite
    cvalsinfinite = [-Inf; Inf]
    # extract other parameters
    nmax          = maxsamplesize(ss, n1)
    p0            = null(params)
    prior         = params.prior
    npriorpivots  = params.npriorpivots
    npivots      = params.npivots
    # get grid for prior and conditional prior
    priorpivots, dcdf   = findgrid(prior, 0, 1, npriorpivots)
    cpriorpivots, dccdf = findgrid(prior, params.pmcrv, 1, 1000) # not performance critical!
    # add type one error rate constraint
    JuMP.@constraint(m,
        sum(dbinom(x1, n1, p0)*_cpr(x1, n1, n, c, p0)*y[x1, n, c] for
            x1 in 0:n1, n in nvals, c in cvals
        ) <= alpha(params)
    )
    # add conditional type two error rate constraint (power)
    for x1 in 0:n1
        posterior1 = dbinom.(x1, n1, cpriorpivots) .* dccdf
        z1 = sum(posterior1)
        JuMP.@constraint(m,
            sum(sum(posterior1 .* _cpr.(x1, n1, n, c, cpriorpivots))/z1*y[x1, n, c] for
                n in nvals, c in cvalsfinite
            ) + sum(y[x1, n, c] for
                n in nvals, c in cvalsinfinite
            ) >= minconditionalpower(params)
        )
        # ensure monotonicity if required
        if (x1 >= 1) & hasmonotoneconditionalpower(params)
            posterior2 = dbinom.(x1 - 1, n1, cpriorpivots) .* dccdf
            z2 = sum(posterior2)
            JuMP.@constraint(m,
                sum(sum(posterior1 .* _cpr.(x1, n1, n, c, cpriorpivots))/z1*y[x1, n, c]
                        - sum(posterior2 .* _cpr.(x1 - 1, n1, n, c, cpriorpivots))/z2*y[x1 - 1, n, c] for
                    n in nvals, c in cvals
                ) >= 0
            )
        end
    end
    # add optimality criterion
    nreq__(p, powerreq) = nreq_(p, powerreq, p0, alpha(params))
    # construct expressions for ROS, upper bound necessary to guarantee finteness!
    JuMP.@expression(m, ros[p in priorpivots],
        1/(params.fs - 1) * sum(
            dbinom(x1, n1, p) * max(0, min(100.0, (n / nreq__(p, 1 - params.beta) - 1))) * y[x1, n, c] for
            x1 in 0:n1, n in nvals, c in cvals
        )
    )
    # construct expressions for power
    JuMP.@expression(m, designpower[p in priorpivots],
        sum(dbinom(x1, n1, p) * _cpr(x1, n1, n, c, p) * y[x1, n, c] for
            x1 in 0:n1, n in nvals, c in cvals
        )
    )
    pivots = [collect(linspace(0, 1 - params.beta, npivots - 1)); 1] # lambda formulaion requires edges!, rup is 0 from powerreq to 1
    JuMP.@variable(m, 0 <= lambda[priorpivots, pivots] <= 1)
    function frup(power, p; rupmax = 100.0)
        if p < params.pmcrv
            return 0.0
        end
        nom   = max(0, nreq__(p, 1 - params.beta) - nreq__(p, power))
        denom = nreq__(p, 1 - params.beta) - nreq__(p, (1 - params.fp)*(1 - params.beta))
        if (denom == 0.0) & (nom > 0)
            return rupmax
        end
        if (denom == 0.0) & (nom == 0)
            return 0.0
        end
        res = max(0, min(rupmax, nom / denom))
        return res
    end
    JuMP.@variable(m, rup[priorpivots])
    for p in priorpivots
        addSOS2(m, [lambda[p, piv] for piv in pivots])
        JuMP.@constraint(m, sum(lambda[p, piv] for piv in pivots) == 1)
        JuMP.@constraint(m, sum(lambda[p, piv]*piv for piv in pivots) == designpower[p]) # defines lambdas!
        JuMP.@constraint(m, sum(lambda[p, piv]*frup(piv, p) for piv in pivots) == rup[p])
    end
    JuMP.@objective(m, Min,
        sum((rup[priorpivots[i]] + ros[priorpivots[i]])*dcdf[i] for i in 1:npriorpivots) # normalized version
    )
    return true
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
    res = min(nmax, max(0, res))
    return res
end
