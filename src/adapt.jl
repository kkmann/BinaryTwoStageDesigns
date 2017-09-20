function adaptstageone(design, data, solver)

  n1           = length(data)
  n1old        = interimsamplesize(design)
  tmp          = data[1:min(n1, n1old)]
  x1star       = tmp[.!DataFrames.isna.(tmp)] |> sum # number of observed resoponses on joint range
  nmissstar    = tmp[DataFrames.isna.(tmp)] |> length # number of NA on joint range
  nmiss        = data[DataFrames.isna.(data[1:n1])] |> length

  params = parameters(design)
  p0     = null(params)

  # reconstruct original problem
  params.samplespace.specialcvalues = [params.samplespace.specialcvalues; n1:(maxsamplesize(params.samplespace) - 1)]
  !possible(n1, samplespace(params)) ? warn("n1 not compatible with sample space") : nothing
  ipm = IPModel(samplespace(params), n1)
  completemodel(ipm, params, n1)

  m = ipm.m
  y = ipm.y

  if (n1 == interimsamplesize(design)) & (nmiss == 0) # nothing to do!
    return(design)
  end
  # fix missingness (locally constant second stage on possible stage one outcomes)
  if nmiss > 0
    for x1_ in (x1star + 1):(x1star + nmiss), n in ipm.nvals, c in ipm.cvals
      if !((c - 1) in ipm.cvals) # c - 1 not available, cannot produce locally constant second stage
        JuMP.@constraint(m,
          y[x1_, n, c] == 0 
        )  
      else # c must be increasing with x1 !
        JuMP.@constraint(m,
          y[x1_, n, c] - y[x1_ - 1, n, c - 1] >= 0 
        )
      end
    end
  end
  # add conditional type one error rate constraint
  if n1 < interimsamplesize(design) # underrunning
    for x1_ in x1star:(x1star + nmissstar)
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
    for x1_ in x1star:(x1star + nmissstar)
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

function adaptstagetwo(design, data, solver)

  n1 = interimsamplesize(design)
  stageonedata = data[1:n1]
  x1 = stageonedata[.!DataFrames.isna.(stageonedata)] |> sum
  nmiss1 = DataFrames.isna.(stageonedata) |> sum
  nnew = length(data)
  nold = samplesize(design, x1)
  for xx1 in x1:(x1 + nmiss1)
    if samplesize(design, xx1) != nold
      error("sample size not constant over possible x1 values, adjust stage one")
    end
  end
  stagetwodata = data[(n1 + 1):end]
  tmp          = data[(n1 + 1):min(nnew, nold)] # joint range
  x2star       = tmp[.!DataFrames.isna.(tmp)] |> sum # number of observed resoponses on joint range
  nmiss2star   = DataFrames.isna.(tmp) |> sum # number of NA on joint range

  println(length(tmp))
  println(x2star)
  println(nmiss2star)

  if nnew == nold # nothing to do!
    return(design)
  end

  params = parameters(design)
  p0     = null(params)
  
  # add observed value of n to nvals
  params.samplespace.specialnvalues = [params.samplespace.specialnvalues; nnew]
  # add all possible c values (required for missing values)
  params.samplespace.specialcvalues = [params.samplespace.specialcvalues; n1:(maxsamplesize(params.samplespace) - 1)]

  # reconstruct original problem
  !possible(n1, samplespace(params)) ? warn("n1 not compatible with sample space") : nothing
  ipm = IPModel(samplespace(params), n1)
  completemodel(ipm, params, n1)

  m = ipm.m
  y = ipm.y

  # constrain n to observed value
  for xx1 in x1:(x1 + nmiss1)
    JuMP.@constraint(m, sum(y[xx1, nnew, c] for c in ipm.cvals) == 1)
  end

  # fix missingness (locally constant second stage on possible stage one outcomes)
  if nmiss1 > 0
    for x1_ in (x1 + 1):(x1 + nmiss1), n in ipm.nvals, c in ipm.cvals
      if !((c - 1) in ipm.cvals) # c - 1 not available, cannot produce locally constant second stage
        JuMP.@constraint(m,
          y[x1_, n, c] == 0 
        )  
      else # c must be increasing with x1 !
        JuMP.@constraint(m,
          y[x1_, n, c] - y[x1_ - 1, n, c - 1] >= 0 
        )
      end
    end
  end

  # add conditional type one error rate constraint 
  if nnew < nold # underrunning
    dn2 = nold - nnew # (positive) difference in sample sizes
    for xx1 in x1:(x1 + nmiss1)
      for x2 in 0:(x2star + nmiss2star)
        # the conditional error of the new design is either 1 or 0
        # but must be smaller or equal than the CE of the old design,
        # therefore the null can only be rejected if the old design would
        # always reject the null as well (CE of old: 1)
        x   = xx1 + x2 .+ 0:dn2 # possible outcomes under old design
        for c in ipm.cvals
          if any(x .>= c) # any non-rejection, c not valid (n must be oberserved value anyway, cf above)
            JuMP.@constraint(m, y[xx1, nnew, c] <= 0)
          end
        end
      end
    end
  end
  if nnew >= nold # overrunning, equal
    dn2 = nnew - nold
    for xx1 in x1:(x1 + nmiss1)
      for x2 in 0:(x2star + nmiss2star)
        # the conditional error of the old design is either 1 or 0
        # and the new design's CE must be smaller or equal.
        # Therefore, no constraint is imposed if the old design already
        # rejected the null (CE of 1) but if it accepted it (CE of 0)
        # the new design must as well
        println()
        println(nnew)
        println(dn2)
        println(xx1 + x2)
        println(criticalvalue(design, xx1))
        if xx1 + x2 <= criticalvalue(design, xx1) # old design does not reject null, new must also (continuation not sensible)
          for c in ipm.cvals
            if xx1 + x2 + dn2 > c # best case: all further samples successes, must still accept null
              JuMP.@constraint(m, y[xx1, nnew, c] <= 0)
            end
          end
        end
      end
    end
  end 

  JuMP.setsolver(ipm.m, solver)
  JuMP.solve(ipm.m) in (:Optimal, :UserLimit) ? nothing : error("no feasible solution reached")
  design = extractsolution(ipm, params)
  _isfeasible(design, params) ? nothing : error("solution infeasible")
  return design
  
end
