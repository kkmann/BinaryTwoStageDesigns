mutable struct MBESS2{TI<:Integer,TR<:Real} <: PointAlternative

  label::String
  samplespace::SampleSpace
  p0::TR
  pmcrv::TR
  prior
  mtoer::TR
  beta::TR
  minconditionalpower::TR
  MONOTONECONDITIONALPOWER::Bool
  npriorpivots::TI # number of pivots for prior evaluation

  function MBESS2{TI,TR}(
    label::String,
    samplespace::SampleSpace,
    p0::TR, pmcrv::TR, prior,
    alpha::TR, beta::TR,
    minconditionalpower::TR, MONOTONECONDITIONALPOWER::Bool,
    npriorpivots::TI
  ) where {TI<:Integer,TR<:Real}

    @checkprob p0 pmcrv alpha beta minconditionalpower
    z = quadgk(prior, 0, 1, abstol = .0001)[1]
    abs(z - 1) > .001 ? 
      error(@sprintf("prior must integrate to one, is %.6f", z)) : nothing
    npriorpivots < 5 ? 
      error(@sprintf("npriorpivots must be at least 5, is %i", npriorpivots)) : nothing
    new(
        label,
        samplespace, 
        p0, pmcrv, prior, 
        alpha, beta, 
        minconditionalpower, MONOTONECONDITIONALPOWER, npriorpivots
    )
    
  end # inner constructor

end # MBESS2


function MBESS2(
  samplespace::SampleSpace,
  p0::Real, pmcrv::Real, prior;
  alpha::Real = 0.05, beta::Real = 0.2,
  minconditionalpower::Real = 0.0, MONOTONECONDITIONALPOWER::Bool = false,
  npriorpivots::Integer = 50,
  label::String = ""
)

  # unify types
  p0, pmcrv,alpha, beta, minconditionalpower = 
    promote(p0, pmcrv, alpha, beta, minconditionalpower)    

  MBESS2{typeof(npriorpivots),typeof(p0)}(
    label,
    samplespace,
    p0, pmcrv, prior,
    alpha, beta,
    minconditionalpower, MONOTONECONDITIONALPOWER,
    npriorpivots
  )

end # external constructor


Base.show(io::IO, params::MBESS2) = print("MBESS2")

maxsamplesize(params::MBESS2) = maxsamplesize(params.samplespace)

isgroupsequential(params::MBESS2) = isgroupsequential(params.ss)

hasmonotoneconditionalpower(params::MBESS2) = params.MONOTONECONDITIONALPOWER

minconditionalpower(params::MBESS2) = params.minconditionalpower

isunimodal(params::MBESS2) = true

prior(params::MBESS2, p) = params.prior(p)




function score(design::Design, params::MBESS2, p::T) where {T<:Real}

    @checkprob p
    ess = SampleSize(design, p) |> mean
    return ess

end # score


function score(design::Design, params::MBESS2)
    return quadgk(
        p -> prior(params, p) * score(design, params, p),
        0, 1, reltol = .001
    )[1]
end


function completemodel(ipm::IPModel, params::MBESS2, n1::Integer)

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
    p0            = params.p0
    pmcrv         = params.pmcrv
    prior         = params.prior
    npriorpivots  = params.npriorpivots
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
        if (x1 >= 1) & hasmonotoneconditionalpower(params)
            JuMP.@constraint(m, # ensure monotonicity if required
                cep[x1] - cep[x1 - 1] >= 0
            )
        end
    end
    z = zeros(n1 + 1) # normalizing constant for prior conditional on X1=x1
    for x1 in 0:n1
        z[x1 + 1] = quadgk(p -> dbinom(x1, n1, p) * prior(p), 0, 1, abstol = .001)[1]
    end
    JuMP.@expression(m, ec[x1 in 0:n1], # expected sample size conditional on X1=x1
        sum(
            n * y[x1, n, c]
            for n in nvals, c in cvals
        )
    )

    JuMP.@constraint(m, # power on conditional prior quantile
        sum(
            dbinom(x1, n1, pmcrv) * _cpr(x1, n1, n, c, pmcrv) * y[x1, n, c]
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
    
    return true

end # completemodel


# ToDo: implement!
_isfeasible(design::Design, params::MBESS2) = true
