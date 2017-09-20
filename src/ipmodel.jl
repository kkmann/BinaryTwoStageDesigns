mutable struct IPModel

  y # the binary assignment variables
  m # the JuMP model
  nvals
  cvalsfinite
  cvals
  n1
  ss

  function IPModel(ss::SampleSpace, n1::Integer)
    
    # plausibility checks
    n1 <= 1           ? error("IPModel: n1 must be >= 1")                    : nothing
    !possible(n1, ss) ? warn("IPModel: n1 not compatible with sample space") : nothing
    
    # extract parameters 
    nmax          = maxsamplesize(ss, n1)
    nvals         = getnvals(ss, n1)
    cvalsfinite   = getcvals(ss, n1)
    cvalsinfinite = [-Inf; Inf]
    cvals         = [-Inf; cvalsfinite; Inf]
    nmincont      = max(n1, ss.nmincont)
    m             = JuMP.Model()

    # indicator variables y[x1, n, c] == 1 iff n(x1) = n, c(x1) = c
    JuMP.@variable(m, y[x1 = 0:n1, n = nvals, c = cvals], Bin)

    for x1 in 0:n1
      # 1) SOS1 constraint (exactly one nonzero entry per x1)
      # JuMP.addSOS1(m, reshape(y[x1, nvals, cvals], length(nvals)*length(cvals)))
      JuMP.@constraint(m,
        sum(y[x1, n, c] for n in nvals, c in cvals) == 1
      )

      # 2) c = -inf (early rejection of null) => n = nmincont (minimal sample size upon continuation)
      # implicitly assumes that smaller sample sizes are better - sensible!
      JuMP.@constraint(m,
        y[x1, nmincont, -Inf] - sum(y[x1, n, -Inf] for n in nvals) >= 0
      )

      # 3) c = inf (early futility stopping) => n = n1 (no second stage)
      # implicitly assumes that smaller sample sizes are better - sensible!
      JuMP.@constraint(m,
        y[x1, n1, Inf] - sum(y[x1, n, Inf] for n in nvals) >= 0
      )

      # 4) contiguous stopping for futility, if c(x1) = inf => c(x1 - 1) = inf
      if x1 > 0
        JuMP.@constraint(m,
          y[x1 - 1, n1, Inf] - y[x1, n1, Inf] >= 0
        )
      end

      # 5) contiguous stopping for efficacy, if c(x1) = -inf => c(x1 + 1) = -inf
      if (x1 < n1)
        JuMP.@constraint(m,
          y[x1 + 1, nmincont, -Inf] - y[x1, nmincont, -Inf] >= 0
        )
      end

      # 6) c(x1) finite => n > n1
      JuMP.@constraint(m,
           (1 - sum(y[x1, n1, c] for c in cvalsfinite)) # 1 iff |c(x1)| == inf and n(x1) > n1
         - sum(y[x1, n, c] for n in nvals, c in cvalsfinite) # 1 iff |c(x1)| < inf
        >= 0
      )

      # 7) c(x1) finite => n - n1 > c - x1 (or |c(x1)| < inf  => n - n1 - c + x1 >= 1)
      JuMP.@constraint(m,
           sum((n - n1 - c + x1) * y[x1, n, c] for n in nvals, c in cvalsfinite)  # n - n1 - c + x1
        >= sum(y[x1, n, c] for n in nvals, c in cvalsfinite) # 1 iff |c(x1)| < inf
      )
      # 7.2) c(x1) finite => c - x1 >= 0 (stage two not already fulfilled) 
      JuMP.@constraint(m,
            sum((c - x1) * y[x1, n, c] for n in nvals, c in cvalsfinite)  # c - x1 if c is finite
        >= (sum(y[x1, n, c] for n in nvals, c in cvalsfinite) - 1) # 0 if |c(x1)| < inf, -1 otherwise
      )

      # 8) n(x1) = n1 or nmincont => |c(x1)| = inf
      JuMP.@constraint(m,
           sum(y[x1, n1, c] + y[x1, nmincont, c] for c in cvalsinfinite)  
        >= sum(y[x1, n1, c] + y[x1, nmincont, c] for c in cvals) # 1 iff n(x1) is n1 or nmincont
      )

    end

    # 9) sample space specific constraints
    JuMP.@constraint(m,
      sum(Real(!possible(n1, n, c, ss)) * y[x1, n, c] for x1 in 0:n1, n in nvals, c in cvalsfinite) == 0 # no impossible configuration should be assigned
    )

    # Group sequential? second stage must be locally constant
    if isgroupsequential(ss)
      for x1 in 1:n1, n in nvals, c in cvalsfinite
        JuMP.@constraint(m, 
          sum(y[x1, nn, c] for nn in nvals, c in cvalsinfinite) + y[x1, n, c] >= y[x1 - 1, n, c] # either early stopping or same n/c
        )
      end
    else # unimodal n
      JuMP.@variable(m, _is_mode[x1 = 0:n1], Bin)
      for x1 in 0:n1
        for x1_ in 1:x1
          JuMP.@constraint(m,
            sum(n*(y[x1_, n, c] - y[x1_ - 1, n, c]) for
              n in nvals[2:end],
              c in cvals
            ) - 3*nmax*_is_mode[x1] >= -3*nmax
          )
        end
        for x1_ in (x1 + 1):n1
          JuMP.@constraint(m,
            sum(n*(y[x1_, n, c] - y[x1_ - 1, n, c]) for
              n in nvals[2:end],
              c in cvals
            ) + 3*nmax*_is_mode[x1] <= 3*nmax
          )
        end
      end
      JuMP.@constraint(m, # at least one mode (can be on boundary as well!)
        sum(_is_mode[x1] for x1 in 0:n1) >= 1 
      )
    end

    new(y, m, nvals, cvalsfinite, cvals, n1, ss)

  end # constructor

end # IPModel

getnvals(ipm::IPModel, n1::Integer) = getnvals(ipm.ss, n1)
getcvals(ipm::IPModel, n1::Integer) = getcvals(ipm.ss, n1)

function extractsolution(ipm::IPModel, params)

    nvec = zeros(Int64, ipm.n1 + 1)
    cvec = zeros(Float64, ipm.n1 + 1) # need float for +/- Inf
    val::Real = 0.0
    current_best_score::Real = 0.0
    for x1 in 0:ipm.n1
        current_best_score = 0.0
        for n in ipm.nvals
            for c in ipm.cvals
                val = JuMP.getvalue(ipm.y[x1, n, c])
                if val > current_best_score # in the end gives the solution where y[x1, n, c] is closest to 1
                    current_best_score = val
                    nvec[x1 + 1] = n
                    cvec[x1 + 1] = c
                end
            end
        end
    end
    # sometimes the probelm is degenerate and stopping regions are not contiguous,
    # manually enforce contiguous stopping:
    inds_futility = find(cvec .== Inf)
    if length(inds_futility) > 0
      cvec[1:inds_futility[1]] = Inf
    end
    inds_efficacy = find(cvec .== -Inf)
    if length(inds_efficacy) > 0
      cvec[inds_efficacy[1]:(ipm.n1 + 1)] = -Inf
    end
    return Design(nvec, cvec, params)

end
