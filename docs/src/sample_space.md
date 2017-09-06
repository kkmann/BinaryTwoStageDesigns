# The Sample Space

Sample space objects can be used to encode the feasible space of binary 
two-stage designs for optimization.

```@docs
SampleSpace

 SampleSpace(
  n1range,
  nmax::TI;
  n2min::TI = 1, 
  maxnfact::TR = Inf, 
  nmincont::TI = 0, 
  maxvariables::TI = 500000, 
  GS::Bool = false
) where {TI<:Integer,TR<:Real}

interimsamplesizerange(ss::SampleSpace)

maxsamplesize(ss::SampleSpace)

maxsamplesize(ss::SampleSpace{TI,TR}, n1::TI2) where {TI<:Integer,TR<:Real,TI2<:Integer}

isgroupsequential(ss::SampleSpace)
```
