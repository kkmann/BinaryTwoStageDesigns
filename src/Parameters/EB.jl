mutable struct EB{TI<:Integer,TR<:Real} <: VagueAlternative

  label::String
  samplespace::SampleSpace
  p0::TR
  pmcrv::TR
  prior
  gamma::TR
  lambda::TR
  a::TR
  b::TR
  targetpower::TR
  k::TR
  mtoer::TR
  minconditionalpower::TR
  MONOTONECONDITIONALPOWER::Bool
  npriorpivots::TI # number of pivots for prior evaluation
  ngpivots::TI # number of pivots for approximation of g

  function EB{TI,TR}(
    label::String,
    samplespace::SampleSpace,
    p0::TR, pmcrv::TR, prior,
    gamma::TR, lambda::TR,
    a::TR, b::TR, targetpower::TR, k::TR,
    alpha::TR,
    minconditionalpower::TR, MONOTONECONDITIONALPOWER::Bool,
    npriorpivots::TI, ngpivots::TI
  ) where {TI<:Integer,TR<:Real}

    @checkprob p0 pmcrv alpha minconditionalpower targetpower
    z = quadgk(prior, 0, 1, abstol = .0001)[1]
    abs(z - 1) > .001 ? 
      error(@sprintf("prior must integrate to one, is %.6f", z)) : nothing
    gamma < 0 ? 
      error(@sprintf("gamma must be positive, is %.6f", gamma)) : nothing
    lambda < 0 ? 
      error(@sprintf("lambda must be positive, is %.6f", gamma)) : nothing
    lambda / gamma > 500  ? 
      warn(@sprintf("lambda/gamma is large (%.3f), potentially infeasible", lambda/gamma)) : nothing
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
      gamma, lambda,
      a, b, targetpower, k,
      alpha,
      minconditionalpower, MONOTONECONDITIONALPOWER,
      npriorpivots, ngpivots
    )

    end # inner constructor

end # EB


function EB( # default values
    samplespace::SampleSpace,
    p0::Real, pmcrv::Real, prior,
    gamma::Real, lambda::Real;
    a::Real = 1, b::Real = 1, targetpower::Real = .8, k::Real = 1,
    alpha::Real = 0.05,
    minconditionalpower::Real = 0.0, MONOTONECONDITIONALPOWER::Bool = true,
    npriorpivots::Integer = 50, ngpivots::Integer = 15,
    label::String = ""
)

  p0, pmcrv, gamma, lambda, alpha, targetpower, a, b, k, minconditionalpower = 
    promote(p0, pmcrv, gamma, lambda, alpha, targetpower, a, b, k, minconditionalpower)    
  npriorpivots, ngpivots = promote(npriorpivots, ngpivots)

  EB{typeof(ngpivots),typeof(p0)}(
    label,
    samplespace,
    p0, pmcrv, prior,
    gamma, lambda,
    a, b, targetpower, k,
    alpha,
    minconditionalpower, MONOTONECONDITIONALPOWER,
    npriorpivots, ngpivots
  )

end # EB


maxsamplesize(params::EB) = maxsamplesize(params.samplespace)

isgroupsequential(params::EB) = isgroupsequential(params.ss)

hasmonotoneconditionalpower(params::EB) = params.MONOTONECONDITIONALPOWER

minconditionalpower(params::EB) = params.minconditionalpower

Base.show(io::IO, params::EB) = print("EB")


function expectedtransformedpower(
  design::Design, params::Union{EB, MBESS}, x1::TI
) where {TI<:Integer}
    
  x1 < 0 ?
    error(@sprintf("x1 must be positive, is %i", x1)) : nothing
  n1       = interimsamplesize(design)
  x1 > n1 ?
    error(@sprintf("x1 must be smaller than n1, is %i", x1)) : nothing
  pmcrv    = mcrv(params)
  phi(p)   = prior(params, p)
  z        = quadgk(p -> phi(p) * Distributions.pdf(Distributions.Binomial(n1, p), x1), pmcrv, 1)[1]
  omega(p) = p < pmcrv ? 0 : phi(p) * Distributions.pdf(Distributions.Binomial(n1, p), x1) / z # conditional prior for stage 2
  return quadgk(p -> g(params, power(design, x1, p)) * omega(p), pmcrv, 1)[1]

end # expectedtransformedpower


function expectedtransformedpower(design::Design, params::Union{EB, MBESS})

  pmcrv  = mcrv(params)
  phi(p) = prior(params, p)
  n1     = interimsamplesize(design)
  z      = quadgk(phi, pmcrv, 1)[1]
  omega(p) = p < pmcrv ? 0 : phi(p) / z # conditional prior given effect
  return quadgk(p -> g(params, power(design, p)) * omega(p), pmcrv, 1)[1]

end # expectedtransformedpower


function expectedcost(design::Design, params::EB, p::T) where {T<:Real}
    
  @checkprob p
  ess = SampleSize(design, p) |> mean
  return params.gamma * ess * (p + params.k * (1 - p))

end # expectedcost


function g(params::EB, power::T) where {T<:Real}

  @checkprob power
  if params.a == params.b == Inf
    if power >= params.targetpower
      return 1.0
    else
      return 0.0
    end
  else
    return Distributions.cdf(Distributions.Beta(params.a, params.b), power)
  end

end # g


function expectedbenefit(design::Design, params::EB, p::T) where {T<:Real}
    
  @checkprob p
  if p < mcrv(params)
    return 0.0
  else
     return params.lambda * g(params, power(design, p))
  end

end # expectedbenefit


function score(design::Design, params::EB, p::T) where {T<:Real}
 
  @checkprob p
  return -(expectedbenefit(design, params, p) - expectedcost(design, params, p))

end # score


function score(design::Design, params::EB)
  
  return quadgk(p -> params.prior(p)*score(design, params, p), 0, 1, reltol = .001)[1]

end # score


function completemodel(ipm::IPModel, params::EB, n1::T) where {T<:Integer}

    ss = samplespace(params)
    !possible(n1, ss) ? warn("completemodel(EB): n1 and sample space incompatible") : nothing
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
    # get grid for conditional prior
    cpriorpivots, dccdf   = findgrid(prior, pmcrv, 1, npriorpivots)
    # add type one error rate constraint
    JuMP.@constraint(m,
        sum(dbinom(x1, n1, p0)*_cpr(x1, n1, n, c, p0)*y[x1, n, c] for
            x1 in 0:n1, n in nvals, c in cvals
        ) <= mtoer(params)
    )
    z = zeros(n1 + 1) # normalizing constants for prior conditional on p >= mcrv, X1=x1
    for x1 in 0:n1
        z[x1 + 1] = quadgk(p -> prior(p)*dbinom(x1, n1, p), pmcrv, 1, abstol = .001)[1]
    end
    JuMP.@expression(m, cep[x1 in 0:n1], # define expressions for conditional power given X1=x1
        sum(quadgk(p -> prior(p)*dbinom(x1, n1, p)*g(params, _cpr.(x1, n1, n, c, p)), pmcrv, 1, abstol = 0.001)[1] / z[x1 + 1] * y[x1, n, c]
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
            (params.gamma*(x1 + params.k*(n1 - x1)) + # stage one cost
             quadgk(p -> params.gamma*(p + params.k*(1 - p))*(n - n1) * dbinom(x1, n1, p)*prior(p)/z[x1 + 1], 0, 1, abstol = .001)[1]  # stage two expected
            ) *
            y[x1, n, c]
            for n in nvals, c in cvals
        )
    )
    # construct expressions for power
    JuMP.@expression(m, designpower[p in cpriorpivots],
        sum(
            dbinom(x1, n1, p)*_cpr(x1, n1, n, c, p) * y[x1, n, c]
            for x1 in 0:n1, n in nvals, c in cvals
        )
    )
    prob_crv = quadgk(p -> prior(p), mcrv(params), 1)[1]
    if params.a == params.b == Inf
        JuMP.@variable(m, exceedstargetpower[1:npriorpivots], Bin)#pivots = [0; params.targetpower - .001; params.targetpower + .001; 1] # lambda formulaion requires edges! model step function as linear approximation
        for i in 1:npriorpivots
            JuMP.@constraint(m,
                designpower[cpriorpivots[i]] - params.targetpower + (1 - exceedstargetpower[i]) >= .00001
            )
        end
        JuMP.@expression(m, obj,
            sum( ec[x1]*z[x1 + 1] for x1 in 0:n1 ) - prob_crv*params.lambda*sum( exceedstargetpower[i]*dccdf[i] for i in 1:npriorpivots ) # negative score!
        )
    else
        lbound = quantile(Distributions.Beta(params.a, params.b), .025) # TODO: implement a = b = inf
        ubound = quantile(Distributions.Beta(params.a, params.b), .975)
        pivots = [0; collect(linspace(lbound, ubound, max(1, ngpivots - 2))); 1] # lambda formulaion requires edges! exploitn fact that function is locally linear!
        JuMP.@variable(m, 0 <= lambdaSOS2[cpriorpivots, pivots] <= 1)
        JuMP.@variable(m, wpwr[cpriorpivots]) # weighted power
        for p in cpriorpivots
            JuMP.addSOS2(m, [lambdaSOS2[p, piv] for piv in pivots])
            JuMP.@constraint(m, sum(lambdaSOS2[p, piv] for piv in pivots) == 1)
            JuMP.@constraint(m, sum(lambdaSOS2[p, piv]*piv for piv in pivots) == designpower[p]) # defines lambdas!
            if p >= params.pmcrv
                JuMP.@constraint(m, sum(lambdaSOS2[p, piv]*g(params, piv) for piv in pivots) == wpwr[p])
            else
                JuMP.@constraint(m, wpwr[p] == 0)
            end
        end
        JuMP.@expression(m, obj,
            sum( ec[x1]*z[x1 + 1] for x1 in 0:n1 ) - prob_crv*params.lambda*sum( wpwr[cpriorpivots[i]]*dccdf[i] for i in 1:npriorpivots ) # negative score!
        )
    end
    JuMP.@objective(m, Min,
        obj
    )
    return true

end # completemodel


# ToDo: implement!
_isfeasible(design::Design, params::EB) = true
