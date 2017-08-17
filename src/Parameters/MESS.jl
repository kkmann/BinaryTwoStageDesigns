"""
    MESS{T_samplespace<:SampleSpace,TR<:Real} <: PointAlternative

This type represents a set of parameters for finding optimal two-stage designs
minimizing the expected sample size on a point in the parameter space subject to
type one and  two error rate constraints.
"""
mutable struct MESS{TR<:Real} <: PointAlternative

  label::String
  samplespace::SampleSpace
  p0::TR
  p1::TR
  mtoer::TR
  beta::TR
  pess::TR
  minconditionalpower::TR
  MONOTONECONDITIONALPOWER::Bool
  minstoppingforfutility::TR

  function MESS{TR}(
    label::String,
    samplespace::SampleSpace,
    p0::TR, p1::TR,
    alpha::TR, beta::TR,
    pess::TR,
    minconditionalpower::TR,
    MONOTONECONDITIONALPOWER::Bool,
    minstoppingforfutility::TR,
  ) where {TR<:Real}

    @checkprob p0 p1 alpha beta pess minconditionalpower minstoppingforfutility
    p0 >= p1 ? 
      error("p0 must be smaller than p1") : nothing
    new(
      label, samplespace, p0, p1, alpha, beta, pess, minconditionalpower, 
      MONOTONECONDITIONALPOWER, minstoppingforfutility
    )

  end # inner constructor

end # MESS

"""
    MESS{T_samplespace<:SampleSpace}(
        samplespace::T_samplespace,
        p0, p1,
        alpha, beta,
        pess;
        minstoppingforfutility::Real   = 0.0,
        minconditionalpower::Real      = 0.0,
        MONOTONECONDITIONALPOWER::Bool = true
    )

Constructs a parameter object of type MESS with the
given values.

# Parameters

| Parameter    | Description |
| -----------: | :---------- |
| samplespacer | a sample space object |
| p0           | upper boundary of the null hypothesis |
| p1           | point alternative to power on |
| alpha        | maximal tolerable type one error rate |
| beta         | maximal tolerable type two error rate on p1 |
| pess         | response rate under which to minimize expected sample size |
| minstoppingforfutility | minimal probability for stopping for futility under p0 |
| minconditionalpower | minimal conditional power upon continuation to stage two |
| MONOTONECONDITIONALPOWER | if true, the conditional power must be monotonously increasing, this constraint is only relevant if nmax is set very restrictively |

# Return Value

An object of type MESS.

# Examples
```julia-repl
julia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)
julia> params = MESS(ss, .2, .4, .05, .2, .4)
```
"""
function MESS(
  samplespace::SampleSpace,
  p0::Real, p1::Real;
  alpha::Real = 0.05, beta::Real = 0.2,
  pess::Real = p1,
  minstoppingforfutility::Real = 0.0,
  minconditionalpower::Real = 0.0,
  MONOTONECONDITIONALPOWER::Bool = true,
  label::String = ""
)

  p0, p1, alpha, beta, pess, minstoppingforfutility, minconditionalpower = 
    promote(p0, p1, alpha, beta, pess, minstoppingforfutility, minconditionalpower)

  MESS{typeof(p0)}(
    label, samplespace, p0, p1, alpha, beta, pess, minconditionalpower, 
    MONOTONECONDITIONALPOWER, minstoppingforfutility
  )

end # MESS


maxsamplesize(params::MESS) = maxsamplesize(params.samplespace)

isgroupsequential(params::MESS) = isgroupsequential(params.ss)

hasmonotoneconditionalpower(params::MESS) = params.MONOTONECONDITIONALPOWER

minconditionalpower(params::MESS) = params.minconditionalpower

Base.show(io::IO, params::MESS) = print("MESS")


function score(design::Design, params::MESS, p::T) where {T<:Real}
    
  @checkprob p
  n = SampleSize(design, p)
  return mean(n)
    
end # score


function score(design::Design, params::MESS)

  n = SampleSize(design, params.pess)
  return mean(n)

end # score


function completemodel(ipm::IPModel, params::MESS, n1::T) where {T<:Integer}

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
  p1            = params.p1
  # add type one error rate constraint
  JuMP.@constraint(m,
      sum(dbinom(x1, n1, p0)*_cpr(x1, n1, n, c, p0)*y[x1, n, c] for
          x1 in 0:n1,
          n  in nvals,
          c  in cvals
      ) <= mtoer(params)
  )
  # add type two error rate constraint (power)
  JuMP.@constraint(m,
      sum(dbinom(x1, n1, p1)*_cpr(x1, n1, n, c, p1)*y[x1, n, c] for
          x1 in 0:n1,
          n  in nvals,
          c  in cvals
      ) >= 1 - params.beta
  )
  # add conditional type two error rate constraint (power)
  for x1 in 0:n1
      # ensure monotonicity if required
      if x1 >= 1 & hasmonotoneconditionalpower(params)
          JuMP.@constraint(m,
              sum(_cpr(x1, n1, n, c, p1)*y[x1, n, c] - _cpr(x1 - 1, n1, n, c, p1)*y[x1 - 1, n, c] for
                  n  in nvals, c in cvals
              ) >= 0
          )
      end
      JuMP.@constraint(m,
          sum(_cpr(x1, n1, n, c, p1)*y[x1, n, c] for
              n  in nvals,
              c  in cvalsfinite
          ) + sum(y[x1, n, c] for
              n  in nvals,
              c  in cvalsinfinite
          ) >= minconditionalpower(params)
      )
  end
  # add constraint for minimal stopping-for-futility probability
  JuMP.@constraint(m,
      sum(dbinom(x1, n1, p0)*y[x1, n1, Inf] for
          x1 in 0:n1
      ) >= params.minstoppingforfutility
  )
  # add optimality criterion
  JuMP.@objective(m, Min,
      sum(dbinom(x1, n1, params.pess)*n*y[x1, n, c] for
          x1 in 0:n1,
          n  in nvals,
          c  in cvals
      )
  )
  return true

end # completemodel


function _isfeasible(design::Design, params::MESS)

    all(power.(design, linspace(0, null(params))) .<= mtoer(params) + .001) ? nothing : throw(InexactError())
    all(power.(design, linspace(alternative(params), 1)) .>= 1 - mtter(params) - .001) ? nothing : throw(InexactError())
    return true

end # _isfeasible
