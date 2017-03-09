abstract SampleSpace

Base.size(::SampleSpace) = ()
Base.getindex(ss::SampleSpace, i) = ss

interimsamplesizerange(ss::SampleSpace) = error("not implemented")
maxsamplesize(ss::SampleSpace) = error("not implemented")
maxsamplesize(ss::SampleSpace, n1) = error("not implemented")
possible{T<:Integer}(n1::T, ss::SampleSpace) = error("not implemented")
possible{T<:Integer, TR<:Real}(n1::T, n::T, c::TR, ss::SampleSpace) = error("not implemented")




type SimpleSampleSpace{T<:Integer} <: SampleSpace
    n1range::Vector{T}
    nmax::T
    n2min
    maxnfact
    nmincont
    function SimpleSampleSpace(n1range, nmax, n2min, maxnfact, nmincont)
        minimum(n1range) < 1 ? throw(InexactError()) : nothing
        maximum(n1range) > nmax ? throw(InexactError()) : nothing
        new(n1range, nmax, n2min, maxnfact, nmincont)
    end
end
SimpleSampleSpace{T<:Integer}(n1range::Vector{T}, nmax::T; n2min = 1, maxnfact = Inf, nmincont = 0) = SimpleSampleSpace{T}(n1range, nmax, n2min, maxnfact, nmincont)
SimpleSampleSpace{T<:Integer}(n1range::UnitRange{T}, nmax::T; n2min = 1, maxnfact = Inf, nmincont = 0) = SimpleSampleSpace{T}(convert(Vector{T}, n1range), nmax, n2min, maxnfact, nmincont)

interimsamplesizerange(ss::SimpleSampleSpace) = ss.n1range
maxsamplesize(ss::SimpleSampleSpace) = ss.nmax
maxsamplesize(ss::SimpleSampleSpace, n1) = convert(Integer, floor(n1*ss.maxnfact))
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
