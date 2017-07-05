"""
    SampleSpace

Abstract type representing generic sample spaces defining the search space for
finding optimal two stage designs.
"""
abstract SampleSpace

Base.size(::SampleSpace) = ()
Base.getindex(ss::SampleSpace, i) = ss

"""
    interimsamplesizerange(ss::SimpleSampleSpace)

Extract the interim sample size range of a `SampleSpace` object.

# Parameters

| Parameter    | Description |
| -----------: | :---------- |
| ss           | SampleSpace object |

# Return Value

The n1range as integer vector.

# Examples
```julia-repl
julia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)
julia> interimsamplesize(ss)
```
"""
interimsamplesizerange(ss::SampleSpace) = error("not implemented")

"""
    maxsamplesize(ss::SampleSpace)

Extract the maximal sample size of a `SampleSpace` object.

# Parameters

| Parameter    | Description |
| -----------: | :---------- |
| ss           | SampleSpace object |

# Return Value

The maximal overall sample size.

# Examples
```julia-repl
julia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)
julia> maxsamplesize(ss)
```
"""
maxsamplesize(ss::SampleSpace) = error("not implemented")

"""
    maxsamplesize(ss::SampleSpace, n1)

Extract the maximal sample size of a `SampleSpace` object given stage-one sample size n1. This might
differ from the maximal sample size due to the factor `maxnfactor`.

# Parameters

| Parameter    | Description |
| -----------: | :---------- |
| ss           | SampleSpace object |
| n1           | interim sample size |

# Return Value

The maximal overall sample size given n1.

# Examples
```julia-repl
julia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)
julia> maxsamplesize(ss, 15)
```
"""
maxsamplesize(ss::SampleSpace, n1) = error("not implemented")

"""
    possible{T<:Integer}(n1::T, ss::SampleSpace)

Returns true if interim sample size n1 is compatible with sample space ss.

# Parameters

| Parameter    | Description |
| -----------: | :---------- |
| n1           | interim sample size |
| ss           | SampleSpace object |

# Return Value

true if n1 is valid for ss, false otherwise

# Examples
```julia-repl
julia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)
julia> possible(15, ss)
```
"""
possible{T<:Integer}(n1::T, ss::SampleSpace) = error("not implemented")

"""
    possible{T<:Integer, TR<:Real}(n1::T, n::T, c::TR, ss::SampleSpace)

TODO
"""
possible{T<:Integer, TR<:Real}(n1::T, n::T, c::TR, ss::SampleSpace) = error("not implemented")

"""
    isgroupsequential(ss::SampleSpace)

Returns true if ss is restricted to group-sequential sample spaces.

# Parameters

| Parameter    | Description |
| -----------: | :---------- |
| ss           | SampleSpace object |

# Return Value

true if ss is group-sequential, false otherwise

# Examples
```julia-repl
julia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)
julia> isgroupsequential(ss)
```
"""
isgroupsequential(ss::SampleSpace) = error("not implemented")

getnvals(ss::SampleSpace, n1) = error("not implemented")
getcvals(ss::SampleSpace, n1) = error("not implemented")

type SimpleSampleSpace{T<:Integer} <: SampleSpace
    n1range::Vector{T}
    nmax::T
    n2min::T
    maxnfact::Real
    nmincont::T
    maxvariables
    GS::Bool # group sequenttial ?
    specialnvalues # required for additional values always included in nvals
    function SimpleSampleSpace(n1range, nmax, n2min, maxnfact, nmincont, maxvariables, GS)
        minimum(n1range) < 1    ? error("minimal n1 must be >= 1")    : nothing
        maximum(n1range) > nmax ? error("maximal n1 must be <= nmax") : nothing
        maxvariables < 100      ? error("maxvariables must be >= 100") : nothing
        maxnfact = (maxnfact == Inf) ? nmax / minimum(n1range) : maxnfact
        new(sort(n1range), nmax, n2min, maxnfact, nmincont, maxvariables, GS, convert(Vector{T}, zeros(0)))
    end
end

SimpleSampleSpace{T<:Integer}(n1range::Vector{T},
    nmax::T; n2min::T = 1, maxnfact::Real = Inf, nmincont::T = 0, maxvariables::T = 500000, GS::Bool = false) = SimpleSampleSpace{T}(n1range, nmax, n2min, maxnfact, nmincont, maxvariables, GS)
"""
    SimpleSampleSpace{T<:Integer}(n1range, nmax::T; n2min::T = 1, maxnfact::Real = Inf, nmincont::T = 0, maxvariables::T = 500000, GS::Bool = false)

Constructs an object of type `SimpleSampleSpace` defining the search space for
finding optimal two stage designs.

# Parameters

| Parameter    | Description |
| -----------: | :---------- |
| n1range      | possible stage-one sample sizes, can be anything which is convertable to an integer vector |
| n2min        | minimal stage-two sample size |
| nmax         | maximal overall sample size (stage one and two combined) |
| maxnfact     | nmax must be smaller than maxnfact*n1 |
| nmincont     | minimal sample size upon rejection of the null hypothesis |
| maxvariables | (approximate) maximal number of binary variables used when finding optimal design; if the sample space is too big this can be used to find approximate solutions to the optimization problem |
| GS           | flag indicating wheather the design should be group-sequential (constant stage two sample size upon continuation) |

# Return Value

An object of type `SimpleSampleSpace` with the respective parameters.

# Examples
```julia-repl
julia> SimpleSampleSpace(10:25, 100, n2min = 5)
```
"""
SimpleSampleSpace{T<:Integer}(n1range, # unspecific, try to convert to integer vector
    nmax::T; n2min::T = 1, maxnfact::Real = Inf, nmincont::T = 0, maxvariables::T = 500000, GS::Bool = false) = SimpleSampleSpace{T}(convert(Vector{T}, n1range), nmax, n2min, maxnfact, nmincont, maxvariables, GS)

interimsamplesizerange(ss::SimpleSampleSpace) = ss.n1range
maxsamplesize(ss::SimpleSampleSpace) = ss.nmax
maxsamplesize(ss::SimpleSampleSpace, n1) = convert(Integer, min(ss.nmax, floor(n1*ss.maxnfact)))
possible{T<:Integer}(n1::T, ss::SimpleSampleSpace) = n1 in ss.n1range
function possible{T<:Integer,TR<:Real}(n1::T, n::T, c::TR, ss::SimpleSampleSpace)
    res = n1 in ss.n1range # well, n1 must be possible
    res = n <= ss.nmax ? res : false # overall sample size constraint must not
                                     # be violated
    res = n1 <= n ? res : false # n1 cannot be larger than n
    res = n <= ss.maxnfact*n1 ? res : false # n cannot be too large in relation
                                            # to n1
    n2  = n - n1
    if (n2 < ss.n2min) & (c != Inf)
        # n2 too small (also affects "early stopping" for efficacy!)
        if (c == -Inf) & (n1 == n) & (n1 >= ss.nmincont) # fine, early stopping for efficacy
            nothing
        else
            res = false
        end
    end
    res = (n  < ss.nmincont) & (c != Inf) ? false : res # only stopping for futility may violate nmincount
    return res
end
isgroupsequential(ss::SimpleSampleSpace) = ss.GS
function getnvals(ss::SimpleSampleSpace, n1)
    nmax    = maxsamplesize(ss, n1)
    frac    = (nmax - n1 + 1) / sqrt(ss.maxvariables/n1)
    if frac > 1
        warn(@sprintf("thinning nvals by factor %.3f", 1/frac))
    else
        frac = 1
    end
    nvals   = convert(Vector{Int64}, round.(collect(n1:frac:nmax)))
    if !(n1 + ss.nmincont in nvals)
        push!(nvals, n1 + ss.nmincont)
    end
    if !((n1 + ss.n2min) in nvals)
        push!(nvals, n1 + ss.n2min)
    end
    if !(maxsamplesize(ss, n1) in nvals)
        push!(nvals, maxsamplesize(ss, n1))
    end
    for snval in ss.specialnvalues # ensure that all special nvals are contained!
        if !(snval in nvals) & (snval >= n1) & (snval <= nmax)
            push!(nvals, snval)
        end
    end
    return convert(Vector{Int64}, sort(nvals))
end
function getcvals(ss::SimpleSampleSpace, n1)
    nmax        = maxsamplesize(ss, n1)
    frac        = nmax / sqrt(ss.maxvariables/n1)
    if frac > 1
        warn(@sprintf("thinning cvals by factor %.3f", 1/frac))
    else
        frac = 1
    end
    cvalsfinite = round.(collect(0:frac:(nmax - 1)))
    return cvalsfinite
end
