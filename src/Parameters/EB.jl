type EB{T_samplespace<:SampleSpace} <: VagueAlternative
    samplespace::T_samplespace
    p0
    pmcrv
    prior
    alpha
    beta
    mincondpower
    gamma
    lambda
    a
    b
    k
    npriorpivots # number of pivots for prior evaluation
    ngpivots # number of pivots for approximation of g
    MONOTONECONDITIONALPOWER
    function EB(
        samplespace, p0, pmcrv, prior, alpha, beta, mincondpower, gamma,
        lambda, a, b, k, npriorpivots, ngpivots, MONOTONECONDITIONALPOWER
    )
        any(!([alpha; p0; pmcrv; mincondpower] .>= 0.0)) ? throw(InexactError()) : nothing
        any(!([alpha; p0; pmcrv; mincondpower] .<= 1.0)) ? throw(InexactError()) : nothing
        abs(quadgk(prior, 0, 1)[1] - 1) <= .001 ? nothing: throw(InexactError())
        new(
            samplespace, p0, pmcrv, prior, alpha, beta, mincondpower, gamma,
            lambda, a, b, k, npriorpivots, ngpivots, MONOTONECONDITIONALPOWER
        )
    end
end
function EB{T_samplespace<:SampleSpace}(
    samplespace::T_samplespace,
    p0, pmcrv, prior, alpha,
    gamma, lambda;
    beta                           = 0.8,
    a::Real                        = 1,
    b::Real                        = 1,
    k::Real                        = 1,
    npriorpivots::Integer          = 33,
    ngpivots::Integer              = 15,
    minconditionalpower::Real      = 0.0,
    MONOTONECONDITIONALPOWER::Bool = true
)
    EB{T_samplespace}(
        samplespace, p0, pmcrv, prior, alpha, beta, minconditionalpower, gamma,
        lambda, a, b, k, npriorpivots, ngpivots, MONOTONECONDITIONALPOWER
    )
end

maxsamplesize(params::EB) = maxsamplesize(params.samplespace)
isgroupsequential{T_samplespace}(params::EB{T_samplespace}) = isgroupsequential(params.ss)
hasmonotoneconditionalpower{T_samplespace}(params::EB{T_samplespace}) = params.MONOTONECONDITIONALPOWER
minconditionalpower(params::EB) = params.mincondpower



function expectedcost(design::AbstractBinaryTwoStageDesign, params::EB, p::Real)
    ess = SampleSize(design, p) |> mean
    return params.gamma*ess * (p + params.k * (1 - p))
end
function g(params::EB, power)
    power < 0 ? error("power cannot be negative") : nothing
    power > 1 ? error("power cannot exceed 1")    : nothing
    return Distributions.cdf(Distributions.Beta(params.a, params.b), power)
end
function expectedbenefit(design::AbstractBinaryTwoStageDesign, params::EB, p::Real)
    if p < mcrv(params)
        return 0.0
    else
        return params.lambda * g(params, power(design, p))
    end
end
score(design::AbstractBinaryTwoStageDesign, params::EB, p::Real) = -(expectedbenefit(design, params, p) - expectedcost(design, params, p))
function score(design::AbstractBinaryTwoStageDesign, params::EB)
    return quadgk(p -> params.prior(p)*score(design, params, p), 0, 1, reltol = .001)[1]
end

function completemodel{T<:Integer}(ipm::IPModel, params::EB, n1::T)
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
    ngpivots      = params.ngpivots
    # get grid for prior and conditional prior
    priorpivots, dcdf   = findgrid(prior, 0, 1, npriorpivots)
    cpriorpivots, dccdf = findgrid(prior, params.pmcrv, 1, 1000) # not performance critical!
    # add type one error rate constraint
    @constraint(m,
        sum(dbinom(x1, n1, p0)*_cpr(x1, n1, n, c, p0)*y[x1, n, c] for
            x1 in 0:n1, n in nvals, c in cvals
        ) <= alpha(params)
    )
    # add conditional type two error rate constraint (power)
    for x1 in 0:n1
        posterior1 = dbinom.(x1, n1, cpriorpivots) .* dccdf # this must be conditional on $\rho>\rho_mcrv
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
            posterior2 = dbinom.(x1 - 1, n1, cpriorpivots) .* dccdf
            z2 = sum(posterior2)
            @constraint(m,
                sum(sum(posterior1 .* _cpr.(x1, n1, n, c, cpriorpivots))/z1*y[x1, n, c]
                        - sum(posterior2 .* _cpr.(x1 - 1, n1, n, c, cpriorpivots))/z2*y[x1 - 1, n, c] for
                    n in nvals, c in cvals
                ) >= 0
            )
        end
    end
    @expression(m, ec[p in priorpivots],
        sum(
            params.gamma * (p + params.k*(1 - p)) * n * dbinom(x1, n1, p) * y[x1, n, c] for
            x1 in 0:n1, n in nvals, c in cvals
        )
    )
    # construct expressions for power
    @expression(m, designpower[p in priorpivots],
        sum(dbinom(x1, n1, p) * _cpr(x1, n1, n, c, p) * y[x1, n, c] for
            x1 in 0:n1, n in nvals, c in cvals
        )
    )
    lbound = quantile(Distributions.Beta(params.a, params.b), .1)
    ubound = quantile(Distributions.Beta(params.a, params.b), .99)
    pivots = [0; collect(linspace(lbound, ubound, max(1, ngpivots - 2))); 1] # lambda formulaion requires edges! exploitn fact that function is locally linear!
    @variable(m, 0 <= lambdaSOS2[priorpivots, pivots] <= 1)
    @variable(m, wpwr[priorpivots]) # weighted power
    for p in priorpivots
        addSOS2(m, [lambdaSOS2[p, piv] for piv in pivots])
        @constraint(m, sum(lambdaSOS2[p, piv] for piv in pivots) == 1)
        @constraint(m, sum(lambdaSOS2[p, piv]*piv for piv in pivots) == designpower[p]) # defines lambdas!
        if p >= params.pmcrv
            @constraint(m, sum(lambdaSOS2[p, piv]*g(params, piv) for piv in pivots) == wpwr[p])
        else
            @constraint(m, wpwr[p] == 0)
        end
    end
    if params.lambda == 0 # constraint optimization
        @constraint(m, sum(wpwr[priorpivots[i]] for i in 1:npriorpivots) >= 1 - params.beta)
        @objective(m, Min,
            sum(ec[priorpivots[i]]*dcdf[i] for i in 1:npriorpivots) # simply minimize expected cost
        )
    else
        @objective(m, Min,
            -sum( (params.lambda*wpwr[priorpivots[i]] - ec[priorpivots[i]])*dcdf[i] for i in 1:npriorpivots) # negative score!
        )
    end
    return true
end

function _isfeasible(design::BinaryTwoStageDesign, params::EB)
    return true
end
