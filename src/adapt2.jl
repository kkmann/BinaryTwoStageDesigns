function adapt_stage_one(design, data, solver; p0 = nothing, p1 = nothing, lambda1 = 1, lambda2 = 1, grid_resolution = 0.05, mincpr = 0.5)
  
  params = parameters(design)

  n1         = length(data)
  data1      = data[1:n1]
  n1old      = interimsamplesize(design)
  n1star     = min(n1, n1old)
  tmp        = data1[1:n1star]
  x1star     = tmp[.!DataFrames.isna.(tmp)] |> sum # number of observed resoponses on joint range
  nna1sstar  = tmp[DataFrames.isna.(tmp)] |> length # number of NA on joint range
  x1         = data[.!DataFrames.isna.(data1)] |> sum # number of observed responses on entire range
  nna1       = data[DataFrames.isna.(data1[1:n1])] |> length
  na_inds    = find(DataFrames.isna.(data))
  nna        = length(na_inds)
  
  if (n1 == n1old) & (nna1 == 0) # nothing to do!
    return(design)
  end
  
  # construct integer programming model
  params.samplespace.specialcvalues = 
    [ params.samplespace.specialcvalues; 
      0:(maxsamplesize(params.samplespace) - 1)
    ]
  ipm = IPModel(samplespace(params), n1)
  m   = ipm.m
  y   = ipm.y

  # overwrite null and alternative if necessary and 
  # construct grid for conditional error control
  p0        = p0 == nothing ? null(params) : p0
  p1        = p1 == nothing ? mcrv(params) : p1
  null_grid = collect(grid_resolution:grid_resolution:p0)
  alt_grid  = collect(p1:grid_resolution:(1 - grid_resolution))

  # create array of possible stage one outcomes (all possible imputations, 2^(# missing))
  stage_one = repmat(transpose(data), 2^nna)
  if nna > 0
    for i in 0:(2^nna - 1)
      stage_one[i + 1, na_inds] = [parse(Int, bin(i, nna)[j]) for j in 1:nna]
    end
  end
  n_imp = size(stage_one, 1)

  # enforce conditional properties for all possible stage-one outcomes
  JuMP.@variable(ipm.m, max_cpr_violation_1 >= 0)
  for i in 1:n_imp
    data_imp = stage_one[i, :] # imputed stage one data
    # stage two must be invariant under imputation
    if i > 1 # i == 1 corresponds to only 0s imputed -> minimal x1
      x1 = Integer(sum(data_imp))
      for n in ipm.nvals, c in ipm.cvals
        if !((c - 1) in ipm.cvals) # c - 1 not available, cannot produce locally constant second stage
          JuMP.@constraint(m,
            y[x1, n, c] == 0 
          )  
        else # c must be increasing with x1 to guarantee invariant stage two!
          JuMP.@constraint(m,
            y[x1, n, c] - y[x1 - 1, n, c - 1] == 0 
          )
        end
      end
    end

    # unconditional type one error rate control, sufficient condition:
    for p in null_grid
      JuMP.@constraint(m, 
        cpr_new(ipm, data_imp, p) <= cpr_old(design, data_imp, p) 
      )
    end

    # make sure that trial continues if possible
    JuMP.@constraint(m, 
      max_cpr_violation_1 >= cpr_old(design, data_imp, p1) - cpr_new(ipm, data_imp, p1) # max_cpr_violation_1 is larger than the biggest undershooting
    )
  end

  p_grid = collect(grid_resolution:grid_resolution:(1 - grid_resolution))

  # ensure monotonicity of cpr in x1
  for x1 in 1:n1, p in [p0; p1] 
    JuMP.@constraint(m,
      sum(_cpr(x1, n1, n, c, p) * y[x1, n, c] - _cpr(x1 - 1, n1, n, c, p) * y[x1 - 1, n, c] for
        n in ipm.nvals, c in ipm.cvals
      ) >= 0
    )
  end

  # ensure new design maintains type one error as well
  JuMP.@constraint(ipm.m,
    pwr_new(ipm, n1, p0) <= power(design, p0)
  )

  # define auxiliary slack variable greater than maximal undershooting of minimal cpr
  mincpr_old = minimum(power.(design, 0:n1old, p1))
  JuMP.@variable(ipm.m, max_cpr_violation_2 >= 0)
  for x1 in 0:n1 # for any possible outcome, max_cpr_violation is larger than max(mincpr, mincpr_old) - cpr
    cvalsfinite   = ipm.cvals[isfinite.(ipm.cvals)]
    cvalsinfinite = ipm.cvals[.!isfinite.(ipm.cvals)]
    JuMP.@constraint(m,
      max_cpr_violation_2 >=
      max(mincpr, mincpr_old) -
      sum(_cpr(x1, n1, n, c, p1) * y[x1, n, c] for
        n  in ipm.nvals,
        c  in cvalsfinite
      ) - 
      sum(y[x1, n, c] for # deactivate constraint for early stopping
        n  in ipm.nvals,
        c  in cvalsinfinite
      )
    )
  end

  # compute expected absolute power difference old/new
  JuMP.@variable(ipm.m, abs_pwr_diff[p in p_grid] >= 0)
  for p in p_grid
    JuMP.@constraint(ipm.m,
      pwr_new(ipm, n1, p) - power(design, p) <= abs_pwr_diff[p]
    )
    JuMP.@constraint(ipm.m,
      -(pwr_new(ipm, n1, p) - power(design, p)) <= abs_pwr_diff[p]
    )
  end
  p_RV = Distributions.Beta(1, 1)
  JuMP.@expression(ipm.m, avg_abs_pwr_diff,
    sum(abs_pwr_diff[p] * Distributions.pdf(p_RV, p) for p in p_grid)
  )

  # objective: expected quadratic difference of sample size under the two designs
  JuMP.@objective(m, :Min, diff_sample_size(design, ipm, n1) + lambda1 * avg_abs_pwr_diff + lambda2 * max_cpr_violation_2 + lambda2 * max_cpr_violation_1)  

  JuMP.setsolver(ipm.m, solver)
  JuMP.solve(ipm.m) in (:Optimal, :UserLimit) ? nothing : error("no feasible solution reached")
  design = extractsolution(ipm, params)
  return design

end


function adapt_stage_two(design, data; p0 = nothing, p1 = nothing)
  
  params   = parameters(design)
  n1       = interimsamplesize(design)
  data1    = data[1:n1]
  x1       = Integer(data1[.!DataFrames.isna.(data1)] |> sum) # number of observed responses stage one
  nna1     = DataFrames.isna.(data1) |> sum
  na1_inds = find(DataFrames.isna.(data1))
  data2    = data[(n1 + 1):end]
  x2       = Integer(data2[.!DataFrames.isna.(data2)] |> sum) # number of observed responses stage one
  nna2     = DataFrames.isna.(data2) |> sum
  na2_inds = find(DataFrames.isna.(data2))
  nna      = nna2 + nna1
  n        = length(data)
  n2       = n - n1
  nold     = samplesize(design, x1)
  n2old    = nold - interimsamplesize(design) # invariant under imputation assumed
  cold     = criticalvalue(design, x1)
  c2old    = cold - x1 # c2old is invariant under missing in stage one

  p0        = p0 == nothing ? null(params) : p0
  p1        = p1 == nothing ? mcrv(params) : p1
  
  if n == nold # nothing to do!
    return(design)
  end

  new_design = deepcopy(design) 
  
  # adjust to new sample size
  new_design.n[(x1:(x1 + nna1)) + 1] = n

  function conditional_probability_to_reject(n2, c2, data2, p) 
    if n2 <= length(data2)
      x2star = sum(data2[1:n2])
      return x2star > c2 ? 1.0 : 0.0
    else
      x2star = sum(data2)
      if x2star > c2
        return 1.0
      else
        dn2 = n2 - length(data2) # remaining patients
        dx2 = c2 - x2star # need dx2 + 1 responses to have x2star + dx2 + 1 > c2
        return sum(Distributions.pdf.(Distributions.Binomial(dn2, p), (dx2 + 1):dn2))
      end
    end
  end
  
  # create array of possible stage one outcomes (all possible imputations, 2^(# missing))
  stage_two = repmat(transpose(data2), 2^nna2)
  if nna2 > 0
    for i in 0:(2^nna2 - 1)
      stage_two[i + 1, na2_inds] = [parse(Int, bin(i, nna2)[j]) for j in 1:nna2]
    end
  end
  n_imp = size(stage_two, 1)

  function possible_c(c)
    c2 = c - x1 
    if c2 < 0
      return false
    end
    res = falses(n_imp)
    for i in 1:n_imp
      data2   = stage_two[i, :]
      cpr_old = conditional_probability_to_reject(n2old, c2old, data2, p0)
      cpr_new = conditional_probability_to_reject(n2, c2, data2, p0)
      if (cpr_new <= cpr_old) & (sum(dbinom.((c2 + 1):n2, n2, p0)) <= sum(dbinom.((c2old + 1):n2old, n2old, p0)))
        res[i] = true
      else
        res[i] = false
      end
    end
    return all(res)
  end

  c = n
  while possible_c(c - 1)
    c -= 1
    if c == 0
      break
    end
  end
  if c == n
    c = Inf
  end

  new_design.c[(x1:(x1 + nna1)) + 1] = c + (0:nna1)

  return(new_design)

end


function diff_sample_size(design, ipm, n1; alpha = 1, beta = 1) 
  # compute expected relative absolute deviation from sample size of original design
  # for old design n() and new n'()
  # default uses uniform beta binomial distributions, alpha and beta of
  # beta mixture can be adjusted though
  dn     = JuMP.AffExpr(0.0)
  n1_old = interimsamplesize(design)
  X1     = Distributions.BetaBinomial(n1_old, alpha, beta)
  X1_bar = Distributions.BetaBinomial(n1, alpha, beta)
  for x1 in 0:n1_old, x1_bar in 0:n1, n in ipm.nvals, c in ipm.cvals
    JuMP.push!(dn, 
      abs(n - samplesize(design, x1)) / samplesize(design, x1) * Distributions.pdf(X1, x1) * Distributions.pdf(X1_bar, x1_bar), 
      ipm.y[x1_bar, n, c]
    )
  end
  return dn
end

function cpr_new(ipm::IPModel, data, p) 
  # conditional probability to reject of new design

  x1 = Integer(sum(data))
  n1 = length(data)
  cpr = JuMP.AffExpr(0.0)
  for n in ipm.nvals, c in ipm.cvals
    JuMP.push!(cpr, 
      _cpr(x1, n1, n, c, p), 
      ipm.y[x1, n, c]
    )
  end
  return cpr
end

function pwr_new(ipm::IPModel, n1, p) 
  # probability to reject of new design

  X1  = Distributions.Binomial(n1, p)
  x1_range = collect(0:n1)
  pwr = JuMP.AffExpr(0.0)
  for x1 in x1_range, n in ipm.nvals, c in ipm.cvals
    JuMP.push!(pwr, 
      Distributions.pdf(X1, x1) * _cpr(x1, n1, n, c, p),
      ipm.y[x1, n, c]
    )
  end
  return pwr
end

function cpr_old(design::Design, data, p) 
  # conditional probability to reject of old design

  n1old = interimsamplesize(design)
  n1new = length(data)
  
  if n1old == n1new 
    x1 = Integer(sum(data))
    return power(design, x1, p)
  end
  if n1old <= n1new
    x1star  = Integer(sum(data[1:n1old])) # number of responses observed in 
                                          # first stage under original design
    if samplesize(design, x1star) == n1old # trial stops early
      return power(design, x1star, p)  
    else # trial continues to stage two
      if n1new >= samplesize(design, x1star) # outcome under original trial
                                              # known
        x = Integer(sum(data[1:samplesize(design, x1star)])) 
        return x > criticalvalue(design, x1star) ? 1.0 : 0.0 
      else # outcome still not necessarily decided
        n2star = samplesize(design, x1star) - n1new # number of samples to 
                                                    # complete original stage two
        xstar   = Integer(sum(data)) # total responses observed so far
        c2star  = criticalvalue(design, x1star) - xstar # c2star + 1 responses
                                                        # needed in remaining 
                                                        # n2star patients
        return sum(dbinom.((c2star + 1):n2star, n2star, p))
      end
    end   
  end
  if n1old > n1new # premature interim analysis
    dn1 = n1old - n1new
    x1star = Integer(sum(data))
    return sum(dbinom.(0:dn1, dn1, p) .* power.(design, x1star + (0:dn1), p))
  end
end
