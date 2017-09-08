function adapt(design, n1, x1star, solver; nmiss = 0)
  
  params = parameters(design)
  p0     = null(params)

  # reconstruct original problem
  !possible(n1, samplespace(params)) ? warn("n1 not compatible with sample space") : nothing
  ipm = IPModel(samplespace(params), n1)
  completemodel(ipm, params, n1)

  m = ipm.m
  y = ipm.y

  # add conditional type one error rate constraint
  if (n1 == interimsamplesize(design)) & (nmiss == 0) # nothing to do!
    return(design)
  end
  if nmiss > 0
    for x1_ in (x1star + 1):(x1star + nmiss)
      for n in ipm.nvals, c in ipm.cvals
        JuMP.@constraint(m,
          y[x1_, n, c] - y[x1_ - 1, n, c] >= 0 # locally constant on potential imputation outcomes
        )
      end
    end
  end
  if n1 < interimsamplesize(design) # underrunning
    for x1_ in x1star:(x1star + nmiss)
      ce_bound = 0.0
      dn1      = interimsamplesize(design) - n1
      for dx1 = 0:dn1
        ce_bound += dbinom(dx1, dn1, p0) * power(design, x1_ + dx1, p0)
      end
      JuMP.@constraint(m, # for every possible x1_ the conditional error must be smaller than under the old design
           sum(_cpr(x1_, n1, n, c, p0) * y[x1_, n, c] for n in ipm.nvals, c in ipm.cvals) 
        <= ce_bound
      )
    end
  end  
  if n1 > interimsamplesize(design) # overrunning
    for x1_ in x1star:(x1star + nmiss)
      dn1    = n1 - interimsamplesize(design)
      ce_old = power(design, x1_, p0)
      JuMP.@constraint(m, # for every possible x1_ the conditional error must be smaller than under the old design
           sum(dbinom(dx1, dn1, p0) * _cpr(x1_ + dx1, n1, n, c, p0) * y[x1_ + dx1, n, c] for n in ipm.nvals, c in ipm.cvals, dx1 in 0:dn1) 
        <= ce_old
      )
    end
  end  
  if n1 == interimsamplesize(design) # only fix missingness
    for x1_ in x1star:(x1star + nmiss)
      ce_old = power(design, x1_, p0)
      JuMP.@constraint(m, # for every possible x1_ the conditional error must be smaller than under the old design
           sum(_cpr(x1_, n1, n, c, p0) * y[x1_, n, c] for n in ipm.nvals, c in ipm.cvals)
        <= ce_old
      )
    end
  end  

  JuMP.setsolver(ipm.m, solver)
  JuMP.solve(ipm.m) in (:Optimal, :UserLimit) ? nothing : error("no feasible solution reached")
  design = extractsolution(ipm, params)
  _isfeasible(design, params) ? nothing : error("solution infeasible")
  return design

end