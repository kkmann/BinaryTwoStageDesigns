struct SampleSize{T<:Real} <: Distributions.DiscreteUnivariateDistribution

  design::Design
  p::T

  function SampleSize{T}(design::Design, p::T) where {T<:Real}
    
    @checkprob p
    new(design, p)

  end # inner constructor

end # SampleSize

SampleSize(design::Design, p::T) where {T<:Real} = SampleSize{T}(design, p)

design(ss::SampleSize) = ss.design

Base.show(io::IO, d::SampleSize) = print("SampleSize")

Base.size(::SampleSize) = ()

Base.getindex(ss::SampleSize, i) = ss

rand(d::SampleSize) = d.design |> interimsamplesize |> n1 -> Distributions.Binomial(n1, d.p) |> rand |> x1 -> samplesize(d.design, x1)


function pdf(d::SampleSize, n::Int)
    
  res = 0.0
  for x1 in 0:interimsamplesize(d.design)
    if samplesize(d.design, x1) == n
      res += Distributions.pdf(Distributions.Binomial(interimsamplesize(d.design), d.p), x1)
    end
  end
  return res

end # pdf 


minimum(d::SampleSize) = interimsamplesize(d.design)

maximum(d::SampleSize) = maximum(samplesize(d.design))


function quantile(d::SampleSize, q::Real)

    nrange = minimum(d):maximum(d)
    res    = 0
    while Distributions.cdf(d, res) <= q
        res += 1
    end
    return res

end # quantile


function mean(d::SampleSize)

    res = 0.0
    for n in minimum(d):maximum(d)
        res += pdf(d, n)*n
    end
    return(res)

end # mean


function var(d::SampleSize)
    
    res = 0.0
    m   = mean(d)
    for n in minimum(d):maximum(d)
        res += pdf(d, n)*(n - m)^2
    end
    return(res)

end # var
