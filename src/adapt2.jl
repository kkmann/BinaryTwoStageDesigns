function adapt_stage_one(design, data, solver)

  n1         = length(data)
  data1      = data[1:n1]
  n1old      = interimsamplesize(design)
  n1star     = min(n1, n1old)
  tmp        = data1[1:n1star]
  x1star     = tmp[.!DataFrames.isna.(tmp)] |> sum # number of observed resoponses on joint range
  nna1sstar  = tmp[DataFrames.isna.(tmp)] |> length # number of NA on joint range
  x1         = data[.!DataFrames.isna.(data1)] |> sum # number of observed responses on joint range
  nna1       = data[DataFrames.isna.(data1[1:n1])] |> length
  
  if (n1 == n1old) & (nna1 == 0) # nothing to do!
    return(design)
  end

  params = parameters(design)
  p0     = null(params)
  params.samplespace.specialcvalues = [params.samplespace.specialcvalues; n1:(maxsamplesize(params.samplespace) - 1)]
  !possible(n1, samplespace(params)) ? warn("n1 not compatible with sample space") : nothing
  ipm = IPModel(samplespace(params), n1)
  m   = ipm.m
  y   = ipm.y

  # implement constraints for recalculation
  
  # 1) fix missingness (locally constant second stage on possible stage one outcomes)
  if nna1 > 0
    for xx1 in (x1star + 1):(x1star + nna1)
      for n in ipm.nvals, c in ipm.cvals
        if !((c - 1) in ipm.cvals) # c - 1 not available, cannot produce locally constant second stage
          JuMP.@constraint(m,
            y[xx1, n, c] == 0 
          )  
        else # c must be increasing with x1 !
          JuMP.@constraint(m,
            y[xx1, n, c] - y[xx1 - 1, n, c - 1] == 0 
          )
        end
      end
    end
  end
  # 2) add conditional type one error rate constraint
  if n1 < n1old # underrunning
    for xx1 in x1star:(x1star + nna1sstar)
      ce_old = 0.0
      dn1    = n1old - n1
      for dx1 = 0:dn1
        ce_old += dbinom(dx1, dn1, p0) * power(design, xx1 + dx1, p0)
      end
      JuMP.@constraint(m, # for every possible xx1 the conditional error must be smaller than under the old design
           sum(_cpr(xx1, n1, n, c, p0) * y[xx1, n, c] for n in ipm.nvals, c in ipm.cvals) 
        <= ce_old
      )
    end
  end  
  if n1 > n1old # overrunning
    for xx1 in x1star:(x1star + nna1sstar)
      dn1    = n1 - n1old
      ce_old = power(design, xx1, p0)
      JuMP.@constraint(m, # for every possible xx1 the conditional error must be smaller than under the old design
           sum(dbinom(dx1, dn1, p0) * _cpr(xx1 + dx1, n1, n, c, p0) * y[xx1 + dx1, n, c] for n in ipm.nvals, c in ipm.cvals, dx1 in 0:dn1) 
        <= ce_old
      )
      # for numerical reasons, a conditional error of 0 might not always lead to
      # early stopping as well
      if ce_old == 0 # manually enforce early stopping if conditional error is 0
        JuMP.@constraint(m,
          sum(1 - y[xx1 + dx1, n1, Inf] for dx1 in 0:(x1 - x1star)) <= 0 # must stop early!
        )
      end
    end
  end  
  if n1 == n1old # only fix missingness
    for x1_ in x1star:(x1star + nminna1sstar)
      ce_old = power(design, x1_, p0)
      JuMP.@constraint(m, # for every possible x1_ the conditional error must be smaller than under the old design
           sum(_cpr(x1_, n1, n, c, p0) * y[x1_, n, c] for n in ipm.nvals, c in ipm.cvals)
        <= ce_old
      )
    end
  end  

  # define expression for quadratic difference in n
  nfunc_old(phat) = samplesize(design, phat*n1old |> floor |> Integer)
  cfunc_old(phat) = criticalvalue(design, phat*n1old |> floor |> Integer)
  lambda = 10^5 * params.samplespace.nmax # required for preventing differences in c with infinite values
  adjust_inf(x) = !isfinite(x) ? lambda : x
  JuMP.@expression(m, ndiff, 0)
  JuMP.@expression(m, cdiff, 0)
  for x1 in 0:n1, n in ipm.nvals, c in ipm.cvals
    push!(ndiff, (n - nfunc_old(x1/n1))^2, y[x1, n, c])
    push!(cdiff, adjust_inf(c - cfunc_old(x1/n1))^2, y[x1, n, c])
  end

  JuMP.@objective(m, :min, cdiff + ndiff)  

  JuMP.setsolver(ipm.m, solver)
  JuMP.solve(ipm.m) in (:Optimal, :UserLimit) ? nothing : error("no feasible solution reached")
  design = extractsolution(ipm, params)
  return design

end
