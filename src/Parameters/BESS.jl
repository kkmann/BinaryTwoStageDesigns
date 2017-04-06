type BESS{T_samplespace<:SampleSpace} <: VagueAlternative
    samplespace::T_samplespace
    p0
    pmcrv
    prior
    a
    b
    k
    alpha
    beta
    minconditionalpower
    MONOTONECONDITIONALPOWER::Bool
    npriorpivots::Integer # number of pivots for prior evaluation
    ngpivots::Integer # number of pivots for approximation of g

    function BESS(
        samplespace,
        p0, pmcrv, prior,
        a, b, k,
        alpha, beta,
        minconditionalpower, MONOTONECONDITIONALPOWER,
        npriorpivots, ngpivots
    )
        @checkprob p0 pmcrv alpha beta minconditionalpower
        z = quadgk(prior, 0, 1, abstol = .0001)[1]
        abs(z - 1)   > .001 ? error(@sprintf("prior must integrate to one, is %.6f", z)) : nothing
        a            <= 0   ? error(@sprintf("a must be positive, is %.3f", a)) : nothing
        b            <= 0   ? error(@sprintf("b must be positive, is %.3f", b)) : nothing
        npriorpivots < 5    ? error(@sprintf("npriorpivots must be at least 5, is %i", npriorpivots)) : nothing
        ngpivots     < 5    ? error(@sprintf("npivots must be at least 5, is %i", ngpivots)) : nothing
        new(
            samplespace,
            p0, pmcrv, prior,
            a, b, k,
            alpha, beta,
            minconditionalpower, MONOTONECONDITIONALPOWER,
            npriorpivots, ngpivots
        )
    end
end
function BESS{T_samplespace<:SampleSpace}( # default values
    samplespace::T_samplespace,
    p0::Real, pmcrv::Real, prior;
    a::Real = 1, b::Real = 1, k::Real = 1,
    alpha::Real = 0.05, beta::Real = 0.2,
    minconditionalpower::Real = 0.0, MONOTONECONDITIONALPOWER::Bool = true,
    npriorpivots::Integer = 50, ngpivots::Integer = 15
)
    BESS{T_samplespace}(
        samplespace,
        p0, pmcrv, prior,
        a, b, k,
        alpha, beta,
        minconditionalpower, MONOTONECONDITIONALPOWER,
        npriorpivots, ngpivots
    )
end

maxsamplesize(params::BESS) = maxsamplesize(params.samplespace)
isgroupsequential{T_samplespace}(params::BESS{T_samplespace}) = isgroupsequential(params.ss)
hasmonotoneconditionalpower{T_samplespace}(params::BESS{T_samplespace}) = params.MONOTONECONDITIONALPOWER
minconditionalpower(params::BESS) = params.minconditionalpower



function g(params::BESS, power)
    @checkprob power
    return Distributions.cdf(Distributions.Beta(params.a, params.b), power)
end
function score(design::AbstractBinaryTwoStageDesign, params::BESS, p::Real)
    @checkprob p
    ess = SampleSize(design, p) |> mean
    return ess * (p + params.k * (1 - p))
end
function score(design::AbstractBinaryTwoStageDesign, params::BESS)
    return quadgk(
        p -> prior(params, p)*score(design, params, p),
        0, 1, reltol = .001
    )[1]
end

function completemodel{T<:Integer}(ipm::IPModel, params::BESS, n1::T)
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
    pmcrv         = mcrv(params)
    prior         = params.prior
    npriorpivots  = params.npriorpivots
    ngpivots      = params.ngpivots
    # add type one error rate constraint
    @constraint(m,
        sum(dbinom(x1, n1, p0)*_cpr(x1, n1, n, c, p0)*y[x1, n, c] for
            x1 in 0:n1, n in nvals, c in cvals
        ) <= alpha(params)
    )
    z = zeros(n1 + 1) # normalizing constants for prior conditional on p >= mcrv, X1=x1
    for x1 in 0:n1
        z[x1 + 1] = quadgk(p -> prior(p)*dbinom(x1, n1, p), pmcrv, 1, abstol = .001)[1]
    end
    @expression(m, cep[x1 in 0:n1], # define expressions for conditional power given X1=x1
        sum(quadgk(p -> prior(p)*dbinom(x1, n1, p)*_cpr.(x1, n1, n, c, p), pmcrv, 1, abstol = 0.001)[1] / z[x1 + 1] * y[x1, n, c]
            for n in nvals, c in cvals
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
    z = zeros(n1 + 1) # normalizing constant for prior conditional on X1=x1
    for x1 in 0:n1
        z[x1 + 1] = quadgk(p -> dbinom(x1, n1, p) * prior(p), 0, 1, abstol = .001)[1]
    end
    @expression(m, ec[x1 in 0:n1], # weighted expected sample size conditional on X1=x1
        sum(
            (x1 + params.k*(n1 - x1) + quadgk(p -> (p + params.k*(1 - p))*(n - n1) * dbinom(x1, n1, p)*prior(p)/z[x1 + 1], 0, 1, abstol = .001)[1]) * y[x1, n, c]
            for n in nvals, c in cvals
        )
    )
    # construct expected weighted power constraint
    cpriorpivots, dccdf = findgrid(prior, params.pmcrv, 1, npriorpivots) # conditional prior pivots
    @expression(m, designpower[p in cpriorpivots], # power on conditional prior pivots
        sum(
            dbinom(x1, n1, p) * _cpr(x1, n1, n, c, p) * y[x1, n, c]
            for x1 in 0:n1, n in nvals, c in cvals
        )
    )
    # pivots for linear approximation of g
    lbound = quantile(Distributions.Beta(params.a, params.b), .025)
    ubound = quantile(Distributions.Beta(params.a, params.b), .975)
    pivots = [0; collect(linspace(lbound, ubound, max(1, ngpivots - 2))); 1] # lambda formulation requires edges! exploit fact that function is locally linear!
    @variable(m, 0 <= lambdaSOS2[cpriorpivots, pivots] <= 1)
    @variable(m, wpwr[cpriorpivots]) # weighted power
    for p in cpriorpivots
        addSOS2(m, [lambdaSOS2[p, piv] for piv in pivots])
        @constraint(m, sum(lambdaSOS2[p, piv] for piv in pivots) == 1)
        @constraint(m, sum(lambdaSOS2[p, piv]*piv for piv in pivots) == designpower[p]) # defines lambdas!
        if p >= params.pmcrv
            @constraint(m, sum(lambdaSOS2[p, piv]*g(params, piv) for piv in pivots) == wpwr[p])
        else
            @constraint(m, wpwr[p] == 0)
        end
    end
    @constraint(m,
        sum( wpwr[cpriorpivots[i]] * dccdf[i] for i in 1:npriorpivots) >= 1 - params.beta
    )
    prior_probs = zeros(n1 + 1)
    for x1 in 0:n1
        prior_probs[x1 + 1] = quadgk(p -> prior(p)*dbinom(x1, n1, p), 0, 1, abstol = .001)[1]
    end
    @objective(m, Min,
        sum( ec[x1]*prior_probs[x1 + 1] for x1 in 0:n1 )
    )
    return true
end

function _isfeasible(design::BinaryTwoStageDesign, params::BESS)
    return true
end
