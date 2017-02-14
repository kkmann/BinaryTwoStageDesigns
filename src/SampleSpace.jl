abstract SampleSpace

Base.size(::SampleSpace) = ()
Base.getindex(ss::SampleSpace, i) = ss

interimsamplesizerange(ss::SampleSpace) = error("not implemented")
maxsamplesize(ss::SampleSpace) = error("not implemented")
maxsamplesize(ss::SampleSpace, n1) = error("not implemented")
possible{T<:Integer}(n1::T, ss::SampleSpace) = error("not implemented")
possible{T<:Integer}(n1::T, n::T, ss::SampleSpace) = error("not implemented")




type SimpleSampleSpace{T<:Integer} <: SampleSpace
    n1range::Vector{T}
    nmax::T
    n2min
    maxnfact
    function SimpleSampleSpace(n1range, nmax, n2min, maxnfact)
        minimum(n1range) < 1 ? throw(InexactError()) : nothing
        maximum(n1range) > nmax ? throw(InexactError()) : nothing
        new(n1range, nmax, n2min, maxnfact)
    end
end
SimpleSampleSpace{T<:Integer}(n1range::Vector{T}, nmax::T; n2min = 1, maxnfact = Inf) = SimpleSampleSpace{T}(n1range, nmax, n2min, maxnfact)
SimpleSampleSpace{T<:Integer}(n1range::UnitRange{T}, nmax::T; n2min = 1, maxnfact = Inf) = SimpleSampleSpace{T}(convert(Vector{T}, n1range), nmax, n2min, maxnfact)

interimsamplesizerange(ss::SimpleSampleSpace) = ss.n1range
maxsamplesize(ss::SimpleSampleSpace) = ss.nmax
maxsamplesize(ss::SimpleSampleSpace, n1) = convert(Integer, floor(n1*ss.maxnfact))
possible{T<:Integer}(n1::T, ss::SimpleSampleSpace) = n1 in ss.n1range
function possible{T<:Integer}(n1::T, n::T, ss::SimpleSampleSpace)
    res = (n1 in ss.n1range) & (n <= ss.nmax) & (n1 <= n)
    res = n - n1 >= ss.n2min ? res : false
    res = n <= ss.maxnfact*n1 ? res : false
    return res
end
