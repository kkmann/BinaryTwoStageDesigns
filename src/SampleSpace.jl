abstract SampleSpace

Base.size(::SampleSpace) = ()
Base.getindex(ss::SampleSpace, i) = ss

interimsamplesizerange(ss::SampleSpace) = error("not implemented")
maxsamplesize(ss::SampleSpace) = error("not implemented")
maxsamplesize(ss::SampleSpace, n1) = error("not implemented")
possible{T<:Integer}(n1::T, ss::SampleSpace) = error("not implemented")
possible{T<:Integer, TR<:Real}(n1::T, n::T, c::TR, ss::SampleSpace) = error("not implemented")
isgroupsequential(ss::SampleSpace) = error("not implemented")
getnvals(ss::SampleSpace, n1) = error("not implemented")
getcvals(ss::SampleSpace) = error("not implemented")



type SimpleSampleSpace{T<:Integer} <: SampleSpace
    n1range::Vector{T}
    nmax::T
    n2min::T
    maxnfact::Real
    nmincont::T
    stepsize::T
    GS::Bool # group sequenttial ?
    function SimpleSampleSpace(n1range, nmax, n2min, maxnfact, nmincont, stepsize, GS)
        minimum(n1range) < 1    ? error("minimal n1 must be >= 1")    : nothing
        maximum(n1range) > nmax ? error("maximal n1 must be <= nmax") : nothing
        maxnfact = maxnfact == Inf ? nmax / minimum(n1range) : maxnfact
        new(n1range, nmax, n2min, maxnfact, nmincont, stepsize, GS)
    end
end
SimpleSampleSpace{T<:Integer}(n1range::Vector{T},
    nmax::T; n2min::T = 1, maxnfact::Real = Inf, nmincont::T = 0, stepsize::T = 1, GS::Bool = false) = SimpleSampleSpace{T}(n1range, nmax, n2min, maxnfact, nmincont, stepsize, GS)
SimpleSampleSpace{T<:Integer}(n1range, # unspecific, try to convert to integer vector
    nmax::T; n2min::T = 1, maxnfact::Real = Inf, nmincont::T = 0, stepsize::T = 1, GS::Bool = false) = SimpleSampleSpace{T}(convert(Vector{T}, n1range), nmax, n2min, maxnfact, nmincont, stepsize, GS)

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
    nvals = collect(n1:ss.stepsize:maxsamplesize(ss, n1))
    if !(n1 + ss.nmincont in nvals)
        push!(nvals, n1 + ss.nmincont)
    end
    if !((n1 + ss.n2min) in nvals)
        push!(nvals, n1 + ss.n2min)
    end
    if !(maxsamplesize(ss, n1) in nvals)
        push!(nvals, maxsamplesize(ss, n1))
    end
    return sort(nvals)
end
function getcvals(ss::SimpleSampleSpace)
    cvalsfinite = collect(0:(maxsamplesize(ss) - 1))
    return cvalsfinite, [-Inf; cvalsfinite; Inf]
end
