"""
    MBESS{TI<:Integer,TR<:Real} <: VagueAlternative

    MBESS(
      samplespace::SampleSpace,
      p0::Real, pmcrv::Real, prior;
      a::Real = 1, b::Real = 1, targetpower::Real = .8, k::Real = 1,
      alpha::Real = 0.05, beta::Real = 0.2,
      minconditionalpower::Real = 0.0, MONOTONECONDITIONALPOWER::Bool = true,
      npriorpivots::Integer = 50, ngpivots::Integer = 15,
      label::String = ""
    )

This type represents a set of parameters for finding optimal two-stage designs
minimizing expected sample size subject to a constraint on expected (transformed)
power.

# Parameters

| Parameter    | Description |
| -----------: | :---------- |
| samplespace  | a sample space object |
| p0           | upper boundary of the null hypothesis |
| pmcrv        | point alternative to power on |
| prior        | prior distribution on true response rate, expected score under the prior conditional on a relevant effect is minimized |
| a            | scaling parameter of nonlinear power transfrom (set to infinity for thresholding) |
| b            | scaling parameter of nonlinear power transfrom |
| targetpower  | target power for thresholding |
| k            | factor for non-responder costs |
| alpha        | maximal type one error rate |
| beta         | expected transformed power must be larger than 1 - beta |
| minconditionalpower | minimal conditional power upon continuation to stage two |
| MONOTONECONDITIONALPOWER | if true, the conditional power must be monotonously increasing, this constraint is only relevant if nmax is set very restrictively |
| npriorpivots | number of pivot points for integration of prior |
| npivots      | number of pivot points piecewise linear approximation of local score |
| label        | descriptive string for parameter set |
"""
mutable struct MBESS{TI<:Integer,TR<:Real} <: VagueAlternative

  label::String
  samplespace::SampleSpace
  p0::TR
  pmcrv::TR
  prior
  a::TR
  b::TR
  targetpower::TR
  k::TR
  mtoer::TR
  beta::TR
  minconditionalpower::TR
  MONOTONECONDITIONALPOWER::Bool
  npriorpivots::TI # number of pivots for prior evaluation
  ngpivots::TI # number of pivots for approximation of g

  function MBESS{TI,TR}(
    label::String,
    samplespace::SampleSpace,
    p0::TR, pmcrv::TR, prior,
    a::TR, b::TR, targetpower::TR, k::TR,
    alpha::TR, beta::TR,
    minconditionalpower::TR, MONOTONECONDITIONALPOWER::Bool,
    npriorpivots::TI, ngpivots::TI
  ) where {TI<:Integer,TR<:Real}

    @checkprob p0 pmcrv alpha beta minconditionalpower targetpower
    z = quadgk(prior, 0, 1, abstol = .0001)[1]
    abs(z - 1) > .001 ? 
      error(@sprintf("prior must integrate to one, is %.6f", z)) : nothing
    a <= 0 ? 
      error(@sprintf("a must be positive, is %.3f", a)) : nothing
    b <= 0 ? 
      error(@sprintf("b must be positive, is %.3f", b)) : nothing
    npriorpivots < 5 ? 
      error(@sprintf("npriorpivots must be at least 5, is %i", npriorpivots)) : nothing
    ngpivots < 5 ? 
      error(@sprintf("npivots must be at least 5, is %i", ngpivots)) : nothing
    new(
        label,
        samplespace, 
        p0, pmcrv, prior, 
        a, b, targetpower, k, 
        alpha, beta, 
        minconditionalpower, MONOTONECONDITIONALPOWER, npriorpivots, ngpivots
    )
    
  end # inner constructor

end # MBESS


function MBESS(
  samplespace::SampleSpace,
  p0::Real, pmcrv::Real, prior;
  a::Real = 1, b::Real = 1, targetpower::Real = .8, k::Real = 1,
  alpha::Real = 0.05, beta::Real = 0.2,
  minconditionalpower::Real = 0.0, MONOTONECONDITIONALPOWER::Bool = true,
  npriorpivots::Integer = 50, ngpivots::Integer = 15,
  label::String = ""
)

  # unify types
  p0, pmcrv, a, b, targetpower, k, alpha, beta, minconditionalpower = 
    promote(p0, pmcrv, a, b, targetpower, k, alpha, beta, minconditionalpower)    
  npriorpivots, ngpivots = promote(npriorpivots, ngpivots)

  MBESS{typeof(ngpivots),typeof(p0)}(
    label,
    samplespace,
    p0, pmcrv, prior,
    a, b, targetpower, k,
    alpha, beta,
    minconditionalpower, MONOTONECONDITIONALPOWER,
    npriorpivots, ngpivots
  )

end # external constructor


Base.show(io::IO, params::MBESS) = print("MBESS")

maxsamplesize(params::MBESS) = maxsamplesize(params.samplespace)

isgroupsequential(params::MBESS) = isgroupsequential(params.ss)

hasmonotoneconditionalpower(params::MBESS) = params.MONOTONECONDITIONALPOWER

minconditionalpower(params::MBESS) = params.minconditionalpower

targetpower(params::MBESS) = params.targetpower


function g(params::MBESS, power::T) where {T<:Real}

    @checkprob power
    if !isfinite(params.a) | !isfinite(params.b)
        # degenerate transform (point measure on target power)
        return power >= params.targetpower ? 1 : 0
    else
        # smooth case
        return Distributions.cdf(Distributions.Beta(params.a, params.b), power)
    end

end # g


function score(design::Design, params::MBESS, p::T) where {T<:Real}

    @checkprob p
    ess = SampleSize(design, p) |> mean
    return ess * (p + params.k * (1 - p))

end # score


function score(design::Design, params::MBESS)
    return quadgk(
        p -> prior(params, p) * score(design, params, p),
        0, 1, reltol = .001
    )[1]
end


function completemodel(ipm::IPModel, params::MBESS, n1::Integer)

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
    JuMP.@constraint(m,
        sum(dbinom(x1, n1, p0) * _cpr(x1, n1, n, c, p0) * y[x1, n, c] for
            x1 in 0:n1, n in nvals, c in cvals
        ) <= mtoer(params)
    )
    z = zeros(n1 + 1) # normalizing constants for prior conditional on p >= mcrv, X1=x1
    for x1 in 0:n1
        z[x1 + 1] = quadgk(p -> prior(p) * dbinom(x1, n1, p), pmcrv, 1, abstol = .001)[1]
    end
    JuMP.@expression(m, cep[x1 in 0:n1], # define expressions for conditional power given X1=x1
        sum(quadgk(p -> prior(p) * dbinom(x1, n1, p) * _cpr.(x1, n1, n, c, p), pmcrv, 1, abstol = 0.001)[1] / z[x1 + 1] * y[x1, n, c]
            for n in nvals, c in cvals
        )
    )
    for x1 in 0:n1
        JuMP.@constraint(m, # add conditional type two error rate constraint (power)
            cep[x1] >= minconditionalpower(params) * (1 - sum(y[x1, n1, c] for c in cvalsinfinite)) # must be conditional on continuation!
        )
        if x1 >= 1 & hasmonotoneconditionalpower(params)
            JuMP.@constraint(m, # ensure monotonicity if required
                cep[x1] - cep[x1 - 1] >= 0
            )
        end
    end
    z = zeros(n1 + 1) # normalizing constant for prior conditional on X1=x1
    for x1 in 0:n1
        z[x1 + 1] = quadgk(p -> dbinom(x1, n1, p) * prior(p), 0, 1, abstol = .001)[1]
    end
    JuMP.@expression(m, ec[x1 in 0:n1], # weighted expected sample size conditional on X1=x1
        sum(
            (x1 + params.k * (n1 - x1) + quadgk(p -> (p + params.k * (1 - p)) * (n - n1) * dbinom(x1, n1, p) * prior(p) / z[x1 + 1], 0, 1, abstol = .001)[1]) * y[x1, n, c]
            for n in nvals, c in cvals
        )
    )
    if params.a == params.b == 1 # special case 1: easier computation!
        z            = quadgk(p -> prior(p), pmcrv, 1, abstol = .0001)[1]
        cprior_probs = zeros(n1 + 1)
        for x1 in 0:n1
            cprior_probs[x1 + 1] = quadgk(p -> prior(p)*dbinom(x1, n1, p) / z, pmcrv, 1, abstol = .0001)[1]
        end
        JuMP.@constraint(m,
            sum( cep[x1]*cprior_probs[x1 + 1] for x1 in 0:n1 ) >= 1 - params.beta
        )
        prior_probs = zeros(n1 + 1)
        for x1 in 0:n1
            prior_probs[x1 + 1] = quadgk(p -> prior(p)*dbinom(x1, n1, p), 0, 1, abstol = .001)[1]
        end
        JuMP.@objective(m, Min,
            sum( ec[x1]*prior_probs[x1 + 1] for x1 in 0:n1 )
        )
    elseif params.a == params.b == Inf # special case 1: easier computation!
        z       = quadgk(p -> prior(p), pmcrv, 1, abstol = .0001)[1]
        ccdf(p) = quadgk(pp -> prior(pp)/z, pmcrv, p, abstol = .001)[1]
        cquant  = Roots.fzero(p -> (1 - ccdf(p)) - params.targetpower, pmcrv, 1) # conditional prior quantile
        JuMP.@constraint(m, # power on conditional prior quantile
            sum(
                dbinom(x1, n1, cquant) * _cpr(x1, n1, n, c, cquant) * y[x1, n, c]
                for x1 in 0:n1, n in nvals, c in cvals
            ) >= 1 - params.beta
        )
        prior_probs = zeros(n1 + 1)
        for x1 in 0:n1
            prior_probs[x1 + 1] = quadgk(p -> prior(p)*dbinom(x1, n1, p), 0, 1, abstol = .001)[1]
        end
        JuMP.@objective(m, Min,
            sum( ec[x1]*prior_probs[x1 + 1] for x1 in 0:n1 )
        )
    else # construct expected weighted power constraint
        cpriorpivots, dccdf = findgrid(prior, params.pmcrv, 1, npriorpivots) # conditional prior pivots
        JuMP.@expression(m, designpower[p in cpriorpivots], # power on conditional prior pivots
            sum(
                dbinom(x1, n1, p) * _cpr(x1, n1, n, c, p) * y[x1, n, c]
                for x1 in 0:n1, n in nvals, c in cvals
            )
        )
        # pivots for linear approximation of g
        lbound = quantile(Distributions.Beta(params.a, params.b), .025)
        ubound = quantile(Distributions.Beta(params.a, params.b), .975)
        pivots = [0; collect(linspace(lbound, ubound, max(1, ngpivots - 2))); 1] # lambda formulation requires edges! exploit fact that function is locally linear!
        JuMP.@variable(m, 0 <= lambdaSOS2[cpriorpivots, pivots] <= 1)
        JuMP.@variable(m, wpwr[cpriorpivots]) # weighted power
        for p in cpriorpivots
            JuMP.addSOS2(m, [lambdaSOS2[p, piv] for piv in pivots])
            JuMP.@constraint(m, sum(lambdaSOS2[p, piv] for piv in pivots) == 1)
            JuMP.@constraint(m, sum(lambdaSOS2[p, piv] * piv for piv in pivots) == designpower[p]) # defines lambdas!
            JuMP.@constraint(m, sum(lambdaSOS2[p, piv] * g(params, piv) for piv in pivots) == wpwr[p])
        end
        JuMP.@expression(m, ewpwr,
            sum( wpwr[cpriorpivots[i]] * dccdf[i] for i in 1:npriorpivots)
        )
        JuMP.@constraint(m,
            ewpwr >= 1 - params.beta
        )
        prior_probs = zeros(n1 + 1)
        for x1 in 0:n1
            prior_probs[x1 + 1] = quadgk(p -> prior(p)*dbinom(x1, n1, p), 0, 1, abstol = .001)[1]
        end
        JuMP.@objective(m, Min,
            sum( ec[x1] * prior_probs[x1 + 1] for x1 in 0:n1 ) # + 2*nmax*absewpwr
        )
    end
    return true

end # completemodel


# ToDo: implement!
_isfeasible(design::Design, params::MBESS) = true
