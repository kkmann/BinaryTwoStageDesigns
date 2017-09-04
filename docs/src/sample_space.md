# The Sample Space

Sample space objects can be used to encode the feasible space of binary 
two-stage designs for optimization.

```@docs
SampleSpace

interimsamplesizerange(ss::SampleSpace)

maxsamplesize(ss::SampleSpace)

maxsamplesize(ss::SampleSpace{TI,TR}, n1::TI2) where {TI<:Integer,TR<:Real,TI2<:Integer}

isgroupsequential(ss::SampleSpace)
```
