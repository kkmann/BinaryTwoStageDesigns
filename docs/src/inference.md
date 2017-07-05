# Inference

```@contents
Pages = ["inference.md"]
Depth = 3
```

## Point estimation

```@docs
BinaryTwoStageDesignEstimator

estimate{T<:Integer}(estimator::BinaryTwoStageDesignEstimator, x1::T, x2::T)

bias{T<:Real}(estimator::BinaryTwoStageDesignEstimator, p::T)

rmse{T<:Real}(estimator::BinaryTwoStageDesignEstimator, p::T)
```

### Maximum likelihood estimator

```@docs
MaximumLikelihoodEstimator
```

### Unbiased estimator

```@docs
RaoBlackwellizedEstimator
```

### Optimal compatible estimator

```@docs
CompatibleEstimator
```

## P values

```@docs
p{T1<:Integer, T2<:Real}(estimator::BinaryTwoStageDesignEstimator, x1::T1, x2::T1, p0::T2)

incompatibleoutcomes(estimator::BinaryTwoStageDesignEstimator)
```

## Confidence intervals

```@docs
ConfidenceInterval

limits

coverage

meanwidth

meaninterval

findinconsistencies
```

### Naive Clopper-Pearson confidence interval

```@docs
NaiveClopperPearsonConfidenceInterval
```

### Clopper-Pearson confidence interval

```@docs
ClopperPearsonConfidenceInterval
```

### Minimum mean width confidence interval

```@docs
MinimumMeanWidthConfidenceInterval
```
