# Binary Two-Stage Designs

```@docs
Design{T1<:Integer, T2<:Real, PType<:Parameters}

Design(n, c, params::PType) where {PType<:Parameters}

Design(n, c)

convert(::Type{DataFrames.DataFrame}, design::Design)

interimsamplesize(design::Design)

parameters(design::Design)

samplesize(design::Design)

samplesize(design::Design, x1::T) where {T<:Integer}

criticalvalue(design::Design)

criticalvalue(design::Design, x1::T) where {T<:Integer}

power{T<:Real}(design::Design, p::T)

power(design::Design, x1::T1, p::T2) where {T1<:Integer, T2<:Real}

expectedpower(design::Design, x1::T, prior::Function; mcrv::Real = mcrv(parameters(design))) where {T<:Integer}

expectedpower(design::Design, prior::Function; mcrv::Real = mcrv(parameters(design)))

stoppingforfutility{T<:Real}(design::Design, p::T)

stoppingforefficacy{T<:Real}(design::Design, p::T)

score(design::Design, params::Parameters)

test(design::Design, x1::T, x2::T) where {T<:Integer}

simulate(design::Design, p::T2, nsim::T1) where {T1<:Integer, T2<:Real}

jeffreysprior(design::Design)

pdf(design::Design, x1::T1, x2::T1, p::T2) where {T1<:Integer, T2<:Real}

save(filename::String, design::Design)

writecsv(filename::String, design::Design; label::String = "") 
```
