"""
    ECPInterval <: ConfidenceInterval

    ECPInterval{T<:Real}(
        estimator::Estimator;
        confidence::T = .9,
        k::Integer = 1001
    )

Exact Clopper-Pearson type confidence interval based on ordering induced by `estimator`.

# Parameters

| Parameter    | Description |
| -----------: | :---------- |
| estimator    | estimator object defining the sample space ordering |
| confidence   | confidence level of the interval |
| k            | number of equally spaced grid-points for invertign the test |

# Examples
```julia-repl
julia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)
julia> interimsamplesize(ss)
julia> design = getoptimaldesign(15, params, Gurobi.GurobiSolver())
julia> est = MaximumLikelihoodEstimator(design, Gurobi.GurobiSolver())
julia> ci = ECPInterval(est, confidence = .9)
```
"""
struct ECPInterval{TE<:Estimator,TR<:Real} <: ConfidenceInterval

  estimator::TE
  confidence::TR
  limits::Array{TR,2}

  function ECPInterval{TE,TR}(
    estimator::TE,
    confidence::TR,
    k::Integer #  = 1001
  ) where {TE<:Estimator,TR<:Real}

    @checkprob confidence
    d      = design(estimator)
    n1     = interimsamplesize(d)
    nmax   = maximum(samplesize(d))
    limits = fill(NaN, n1 + 1, nmax - n1 + 1, 2)
    supp   = support(d)

    n1     = interimsamplesize(d)
    nmax   = maximum(samplesize(d))
    grid   = collect(linspace(0, 1, k))
    supp   = support(d) # compute p values by hand is faster than calling p()
    tmp    = zeros(Float64, size(supp, 1), k)
    estimates = estimate.(estimator, supp[:, 1], supp[:, 2])
    for i in 1:size(supp, 1)
        tmp[i, :] = pdf.(d, supp[i, 1], supp[i, 2], linspace(0, 1, k))
    end
    probs = [supp samplesize.(d, supp[:, 1]) estimates tmp]
    probs = sortrows(probs, by = x -> x[4], rev = true) # sort in descending order of estimates
    # compute p vals for H0: p < p0
    pvals_leq = probs |> x -> cumsum(x, 1)
    pvals_leq[:, 1:4] = probs[:, 1:4]

    # compute pvals for H0: p > p0
    pvals_geq = sortrows(probs, by = x -> x[4]) |> x -> cumsum(x, 1)
    pvals_geq[:, 1:4] = sortrows(probs, by = x -> x[4])[:, 1:4]
    pvals_geq = sortrows(pvals_geq, by = x -> x[4], rev = true)

    for i in 1:size(supp, 1)
        x1, x2 = supp[i, :]
        limits[x1 + 1, x2 + 1, :] = get_limits(d, x1, x2, confidence, pvals_leq, pvals_geq, grid, k)
    end
    new(estimator, confidence, limits)

  end # inner constructor

end # ECPInterval


function ECPInterval(
  estimator::TE; 
  confidence::TR = .9, k::Integer = 1001
) where {TE<:Estimator,TR<:Real} 

  return ECPInterval{TE,TR}(estimator; confidence, k)

end


estimator(ci::ECPInterval) = ci.estimator

design(ci::ECPInterval) = ci |> estimator |> design


function get_limits(
  d::Design, x1::T1, x2::T1, confidence::T2, pvals_leq, pvals_geq, grid, k::T1
) where {T1<:Integer, T2<:Real}

  pvals_leq = pvals_leq[(pvals_leq[:, 1] .== x1) & (pvals_leq[:, 2] .== x2), 5:(k + 4)][1, :]
  pvals_geq = pvals_geq[(pvals_geq[:, 1] .== x1) & (pvals_geq[:, 2] .== x2), 5:(k + 4)][1, :]

  indlow = findfirst(pvals_leq .>= (1 - confidence)/2.0)
  indlow = indlow == 0 ? 1 : indlow
  indup  = findlast(pvals_geq .>= (1 - confidence)/2.0)
  indup  = indup == 0 ? k : indup
  limits = [grid[indlow]; grid[indup]]
  if indlow == 1
      limits[1] = 0.0
  end
  if indup == k
      limits[2] = 1.0
  end
  if x1 + x2 == 0
      limits[1] = 0.0
  end
  if x1 + x2 == samplesize(d, x1)
      limits[2] = 1.0
  end
  return limits

end


function limits(ci::ECPInterval, x1::T, x2::T) where {T<:Integer}

  ispossible(design(ci), x1, x2) ? nothing : error("x1 x2 combination not possible")
  return convert(Vector{Float64}, view([ci.limits[x1 + 1, x2 + 1, 1] ci.limits[x1 + 1, x2 + 1, 2]], 1:2))

end
