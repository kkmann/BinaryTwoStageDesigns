# Binary Two-Stage Designs

```@docs
BinaryTwoStageDesign{T1<:Integer, T2<:Real, PType<:Parameters}

convert(::Type{DataFrames.DataFrame}, design::BinaryTwoStageDesign)

interimsamplesize(design::BinaryTwoStageDesign)

parameters(design::BinaryTwoStageDesign)

samplesize(design::BinaryTwoStageDesign)

criticalvalue(design::BinaryTwoStageDesign)

power{T<:Real}(design::BinaryTwoStageDesign, p::T)

stoppingforfutility{T<:Real}(design::BinaryTwoStageDesign, p::T)

test{T<:Integer}(design::BinaryTwoStageDesign, x1::T, x2::T)

simulate{T1<:Integer, T2<:Real}(design::BinaryTwoStageDesign, p::T2, nsim::T1)

score(design::BinaryTwoStageDesign, params::Parameters)

jeffreysprior(design::BinaryTwoStageDesign)

pdf{T1<:Integer, T2<:Real}(design::BinaryTwoStageDesign, x1::T1, x2::T1, p::T2)
```
