# Inference

```@contents
Pages = ["inference.md"]
Depth = 3
```

## Point estimation

```@docs
Estimator

estimate{T<:Integer}(estimator::Estimator, x1::T, x2::T)

bias{T<:Real}(estimator::Estimator, p::T)

rmse{T<:Real}(estimator::Estimator, p::T)
```

### Maximum likelihood estimator

```@docs
MLEstimator(design::Design)
```

### Unbiased estimator

```@docs
RBEstimator(design::Design)
```

### Optimal compatible estimator

```@docs
OCEstimator{TS<:MathProgBase.AbstractMathProgSolver}(
    design::Design,
    solver::TS;
    prior::Function = jeffreysprior(design),
    k = 100
)
```

## P values

```@docs
pvalue(estimator::Estimator, x1::T1, x2::T1, p0::T2) where {T1<:Integer, T2<:Real}

incompatibleoutcomes(estimator::Estimator)
```

## Confidence intervals

```@docs
ConfidenceInterval

limits(ci::ConfidenceInterval, x1::T, x2::T) where {T<:Integer}

coverage(ci::ConfidenceInterval, p::T; orientation::String = "overall") where {T<:Real}

meanwidth(ci::ConfidenceInterval, p::T) where {T<:Real}

meaninterval(ci::ConfidenceInterval, p::T) where {T<:Real}

findinconsistencies(ci::ConfidenceInterval, p0::T) where {T<:Real}
```

### Naive Clopper-Pearson confidence interval

```@docs
CPInterval{T<:Real}(
  design::Design;
  confidence::T = .9
)
```

### Clopper-Pearson confidence interval

```@docs
ECPInterval{T<:Real}(
  estimator::Estimator;
  confidence::T = .9,
  k::Integer = 1001
)
```

### Minimum mean width confidence interval

```@docs
MMWInterval(
  estimator::TE,
  rho0::TR,
  prior::Function,
  solver::MathProgBase.AbstractMathProgSolver;
  confidence::TR = 0.9,
  ngrid::Integer = 100
) where {TE<:Estimator,TR<:Real}
```
