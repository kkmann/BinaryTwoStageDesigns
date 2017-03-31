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
    beta
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
        beta,
        mincondpower,
        gamma, lambda,
        a, b, k,
        smoothness
    )
        any(!([alpha; p0; pmcrv; mincondpower] .>= 0.0)) ? throw(InexactError()) : nothing
        any(!([alpha; p0; pmcrv; mincondpower] .<= 1.0)) ? throw(InexactError()) : nothing
        abs(quadgk(prior, 0, 1)[1] - 1) <= .001 ? nothing: throw(InexactError())
        new(samplespace, p0, pmcrv, prior, alpha, beta, mincondpower, gamma, lambda, a, b, k, smoothness)
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

function expectedcost(design::AbstractBinaryTwoStageDesign, params::EB, p::Real)
    ess = SampleSize(design, p) |> mean
    return params.gamma*ess * (p + params.k * (1 - p))
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


function _createProblem{T<:Integer}(
    n1::T,      # stage one sample size
    params::EB;
    npivots      = 15,
    npriorpivots = 33
)
    ss = samplespace(params)
    nmax = maxsamplesize(ss, n1)
    possible(n1, ss) ? nothing : throw(InexactError())
    a  = alpha(params)
    p0 = null(params)
    prior = params.prior
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
    # priorpivots  = collect(linspace(0.0, 1.0, npriorpivots + 2))[2:(npriorpivots + 1)] # leave out boundary values!
    # dp           = priorpivots[2] - priorpivots[1]
    # priorvals    = prior.(priorpivots) ./ sum(prior.(priorpivots) .* dp) # normalize to 1
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
    pivots = [0; collect(linspace(lbound, ubound, max(1, npivots - 2))); 1] # lambda formulaion requires edges! exploitn fact that function is locally linear!
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
    return m, y
end

function _isfeasible(design::BinaryTwoStageDesign, params::EB)
    return true
end

# utility
function g(params::EB, power)
    power < 0 ? throw(InexactError()) : nothing
    power > 1 ? throw(InexactError()) : nothing
    return Distributions.cdf(Distributions.Beta(params.a, params.b), power)
end
