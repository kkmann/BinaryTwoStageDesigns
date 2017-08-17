"""
    Estimator

Abstract base type for all estimators.
"""
abstract type Estimator end

# these two enable array broadcasting of all methods!
Base.size(::Estimator) = ()

Base.getindex(estimator::Estimator, i) = estimator


design(estimator::Estimator) = try return estimator.design catch error("not implemented") end

"""
    estimate{T<:Integer}(estimator::Estimator, x1::T, x2::T)

Estimate the response rate from observed x1 and x2.

# Parameters

| Parameter    | Description |
| -----------: | :---------- |
| estimator    | any BinaryTwoStageDeisgnEstimator object |
| x1           | stage-one responses |
| x2           | stage-two responses |

# Return Value

Real, estimated response rate.

# Examples
```julia-repl
julia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)
julia> interimsamplesize(ss)
julia> design = getoptimaldesign(15, params, solver = Gurobi.GurobiSolver())
julia> est = MLE(design)
julia> estimate(est, 0, 0)
```
"""
estimate(estimator::Estimator, x1::T, x2::T) where {T<:Integer} = error("not implemented")

"""
    p{T1<:Integer, T2<:Real}(estimator::Estimator, x1::T1, x2::T1, p0::T2)

Compute the p value after observing (x1, x2) for null hypothesis H0: p <= p0 with
respect to ordering induced by `estimator`.

# Parameters

| Parameter    | Description |
| -----------: | :---------- |
| estimator    | any BinaryTwoStageDeisgnEstimator object |
| x1           | stage-one responses |
| x2           | stage-two responses |
| p0           | upper boundary of null hypothesis |

# Return Value

Real, p value.

# Examples
```julia-repl
julia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)
julia> interimsamplesize(ss)
julia> design = getoptimaldesign(15, params, solver = Gurobi.GurobiSolver())
julia> est = MLE(design)
julia> p(est, 0, 0, .2)
```
"""
function pvalue(
  estimator::Estimator, x1::T1, x2::T1, p0::T2
) where {T1<:Integer, T2<:Real}

  checkx1x2(x1, x2, design(estimator))
  supp          = support(design(estimator))
  estimates     = estimate.(estimator, supp[:, 1], supp[:, 2])
  estimateobs   = estimate(estimator, x1, x2)
  ismoreextreme = estimates .>= estimateobs
  return pdf.(design(estimator), supp[ismoreextreme, 1], supp[ismoreextreme, 2], p0) |> sum

end

"""
    bias{T<:Real}(estimator::Estimator, p::T)

Bias of `estimator` given response rate `p`.

# Parameters

| Parameter    | Description |
| -----------: | :---------- |
| estimator    | any BinaryTwoStageDeisgnEstimator object |
| p0           | upper boundary of null hypothesis |

# Return Value

Real, bias given `p`.

# Examples
```julia-repl
julia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)
julia> interimsamplesize(ss)
julia> design = getoptimaldesign(15, params, solver = Gurobi.GurobiSolver())
julia> est = MLE(design)
julia> bias(est, .3)
```
"""
function bias(estimator::Estimator, p::T) where {T<:Real}

  supp      = support(design(estimator))
  estimates = estimate.(estimator, supp[:, 1], supp[:, 2])
  probs     = pdf.(design(estimator), supp[:, 1], supp[:, 2], p)
  return dot(probs, estimates) - p

end

"""
    rmse{T<:Real}(estimator::Estimator, p::T)

Root mean squared error of `estimator` given response rate `p`.

# Parameters

| Parameter    | Description |
| -----------: | :---------- |
| estimator    | any BinaryTwoStageDeisgnEstimator object |
| p0           | upper boundary of null hypothesis |

# Return Value

Real, RMSE given `p`.

# Examples
```julia-repl
julia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)
julia> interimsamplesize(ss)
julia> design = getoptimaldesign(15, params, solver = Gurobi.GurobiSolver())
julia> est = MLE(design)
julia> rmse(est, .3)
```
"""
function rmse(estimator::Estimator, p::T) where {T<:Real}

  supp     = support(design(estimator))
  return sqrt( sum(
    pdf.(design(estimator), supp[:, 1], supp[:, 2], p) .*
      (estimate.(estimator, supp[:, 1], supp[:, 2]) .- p).^2
  ) )

end

"""
    incompatibleoutcomes(estimator::Estimator)

Find outcomes where the induced p value implies different decision than the
underlying design

# Parameters

| Parameter    | Description |
| -----------: | :---------- |
| estimator    | any BinaryTwoStageDeisgnEstimator object |

# Return Value

Array with respective outcomes

# Examples
```julia-repl
julia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)
julia> interimsamplesize(ss)
julia> design = getoptimaldesign(15, params, solver = Gurobi.GurobiSolver())
julia> est = MLE(design)
julia> incompatibleoutcomes(est)
```
"""
function incompatibleoutcomes(estimator::Estimator)

    supp = support(design(estimator))
    res  = []
    p0   = estimator |> design |> parameters |> null
    a    = estimator |> design |> parameters |> alpha
    for i in 1:size(supp, 1)
        x1, x2 = supp[i, :]
        if (pvalue(estimator, x1, x2, p0) <= a) != test(design(estimator), x1, x2)
            push!(res, supp[i, :])
        end
    end
    return hcat(res...)'
    
end
