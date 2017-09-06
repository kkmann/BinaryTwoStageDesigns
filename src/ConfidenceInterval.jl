"""
    ConfidenceInterval

Abstract base type for all confidence interval types.
"""
abstract type ConfidenceInterval end

# these two enable array broadcasting of all methods!
Base.size(::ConfidenceInterval) = ()

Base.getindex(ci::ConfidenceInterval, i) = ci

"""
    limits{T<:Integer}(ci::ConfidenceInterval, x1::T, x2::T)

Return the confidence interval limits for observed x1 and x2.

# Parameters

| Parameter    | Description |
| -----------: | :---------- |
| ci           | any ConfidenceInterval object |
| x1           | stage-one responses |
| x2           | stage-two responses |

# Return Value

Two element Real vector of limits [lower, upper].

# Examples
```julia-repl
julia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)
julia> interimsamplesize(ss)
julia> design = getoptimaldesign(15, params, solver = Gurobi.GurobiSolver())
julia> est = MLE(design)
julia> ci = ClopperPearsonConfidenceInterval(est, confidence = .9)
julia> limits(ci, 0, 0)
```
"""
limits(ci::ConfidenceInterval, x1::T, x2::T) where {T<:Integer} = error("not implemented!")

confidence(ci::ConfidenceInterval) = try return ci.confidence catch error("not implemented") end

design(ci::ConfidenceInterval) = try return ci.design catch error("not implemented") end

"""
    coverage(ci::ConfidenceInterval, p::T; orientation::String = "overall") where {T<:Real}

Return coverage of given confidence interval and response rate `p`.

# Parameters

| Parameter    | Description |
| -----------: | :---------- |
| ci           | any ConfidenceInterval object |
| p            | response rate |
| orientation  | string indicating the coverage type - "overall", "lower", or "upper" |

# Return Value

Coverage probability given `p`.

# Examples
```julia-repl
julia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)
julia> interimsamplesize(ss)
julia> design = getoptimaldesign(15, params, solver = Gurobi.GurobiSolver())
julia> est = MLE(design)
julia> ci = ClopperPearsonConfidenceInterval(est, confidence = .9)
julia> coverage(ci, .5)
```
"""
function coverage(ci::ConfidenceInterval, p::T; orientation::String = "overall") where {T<:Real} # or upper or lower

  @checkprob p 
  supp   = support(design(ci))
  res    = 0.0
  for i in 1:size(supp, 1)
    x1, x2 = supp[i, :]
    lim = limits(ci, x1, x2)
    if (orientation == "overall") & (lim[1] <= p) & (lim[2] >= p)
        res += pdf(design(ci), x1, x2, p)
    elseif (orientation == "upper") & (lim[2] >= p)
        res += pdf(design(ci), x1, x2, p)
    elseif (orientation == "lower") & (lim[1] <= p)
        res += pdf(design(ci), x1, x2, p)
    end
  end
  return res

end


"""
    meanwidth(ci::ConfidenceInterval, p::T) where {T<:Real}

Mean width of given confidence interval and response rate `p`.

# Parameters

| Parameter    | Description |
| -----------: | :---------- |
| ci           | any ConfidenceInterval object |
| p            | response rate |

# Return Value

Mean width given `p`.

# Examples
```julia-repl
julia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)
julia> interimsamplesize(ss)
julia> design = getoptimaldesign(15, params, solver = Gurobi.GurobiSolver())
julia> est = MLE(design)
julia> ci = ClopperPearsonConfidenceInterval(est, confidence = .9)
julia> meanwidth(ci, .5)
```
"""
function meanwidth(ci::ConfidenceInterval, p::T) where {T<:Real}

  d = design(ci)
  supp   = support(d)
  res    = 0.0
  for i in 1:size(supp, 1)
    x1, x2 = supp[i, :]
    res = res + pdf(d, x1, x2, p)*(limits(ci, x1, x2)[2] - limits(ci, x1, x2)[1])
  end
  return res

end


"""
    meaninterval(ci::ConfidenceInterval, p::T) where {T<:Real}

Mean interval (average limits) of given confidence interval and response rate `p`.

# Parameters

| Parameter    | Description |
| -----------: | :---------- |
| ci           | any ConfidenceInterval object |
| p            | response rate |

# Return Value

Mean interval given `p`.

# Examples
```julia-repl
julia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)
julia> interimsamplesize(ss)
julia> design = getoptimaldesign(15, params, solver = Gurobi.GurobiSolver())
julia> est = MLE(design)
julia> ci = ClopperPearsonConfidenceInterval(est, confidence = .9)
julia> meaninterval(ci, .5)
```
"""
function meaninterval(ci::ConfidenceInterval, p::T) where {T<:Real}

  d = design(ci)
  supp   = support(d)
  mean_limits = [0.0; 0.0]
  for i in 1:size(supp, 1)
      x1, x2      = supp[i, :]
      mean_limits = mean_limits + pdf(d, x1, x2, p)*limits(ci, x1, x2)
  end
  return mean_limits

end


"""
    findinconsistencies(ci::ConfidenceInterval, p0::T) where {T<:Real}

Return outcomes where the confidence interval contradicts the designs test decision.

# Parameters

| Parameter    | Description |
| -----------: | :---------- |
| ci           | any ConfidenceInterval object |
| p            | response rate |

# Return Value

Mean interval given `p`.

# Examples
```julia-repl
julia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)
julia> interimsamplesize(ss)
julia> design = getoptimaldesign(15, params, solver = Gurobi.GurobiSolver())
julia> est = MLE(design)
julia> ci = ClopperPearsonConfidenceInterval(est, confidence = .9)
julia> findinconsistencies(ci)
```
"""
function findinconsistencies(ci::ConfidenceInterval, p0::T) where {T<:Real}

  d = design(ci)
  supp   = support(d)
  res    = []
  for i in 1:1:size(supp, 1)
      x1, x2 = supp[i, :]
      if (limits(ci, x1, x2)[1] <= p0) & (x1 + x2 > d.c[x1 + 1])
          push!(res, [x1 x2 d.c[x1 + 1] limits(ci, x1, x2)[1] p0])
      end
      if (limits(ci, x1, x2)[1] > p0) & (x1 + x2 <= d.c[x1 + 1])
          push!(res, [x1 x2 d.c[x1 + 1] limits(ci, x1, x2)[1] p0])
      end
  end
  res = vcat(res...)
  return res
  
end
