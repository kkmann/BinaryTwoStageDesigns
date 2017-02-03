abstract SampleSpace

Base.size(::SampleSpace) = ()
Base.getindex(ss::SampleSpace, i) = ss

interimsamplesizerange(ss::SampleSpace) = error("not implemented")
maxsamplesize(ss::SampleSpace) = error("not implemented")
possible{T<:Integer}(n1::T, ss::SampleSpace) = error("not implemented")
possible{T<:Integer}(n1::T, n::T, ss::SampleSpace) = error("not implemented")




type SimpleSampleSpace{T<:Integer} <: SampleSpace
    n1range::Vector{T}
    nmax::T
    function SimpleSampleSpace(n1range, nmax)
        minimum(n1range) < 1 ? throw(InexactError()) : nothing
        maximum(n1range) > nmax ? throw(InexactError()) : nothing
        new(n1range, nmax)
    end
end
SimpleSampleSpace{T<:Integer}(n1range::Vector{T}, nmax::T) = SimpleSampleSpace{T}(n1range, nmax)
SimpleSampleSpace{T<:Integer}(n1range::UnitRange{T}, nmax::T) = SimpleSampleSpace{T}(convert(Vector{T}, n1range), nmax)

interimsamplesizerange(ss::SimpleSampleSpace) = ss.n1range
maxsamplesize(ss::SimpleSampleSpace) = ss.nmax
possible{T<:Integer}(n1::T, ss::SimpleSampleSpace) = n1 in ss.n1range
possible{T<:Integer}(n1::T, n::T, ss::SimpleSampleSpace) = (n1 in ss.n1range) & (n <= ss.nmax) & (n1 <= n)
