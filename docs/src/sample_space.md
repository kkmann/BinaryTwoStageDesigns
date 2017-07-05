# Sample Spaces

Sample spaces are used to encode the feasible space of binary two-stage designs
for optimization.
Currently, only a single option `SimpleSampleSpace` exists which allows to
specify a range of possible stage-one sample sizes, maximal overall sample size
and various nicety constraints.

```@docs
SampleSpace

SimpleSampleSpace{T<:Integer}(n1range, nmax::T; n2min::T = 1, maxnfact::Real = Inf, nmincont::T = 0, maxvariables::T = 500000, GS::Bool = false)

```
