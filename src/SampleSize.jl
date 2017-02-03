type SampleSize <: Distributions.DiscreteUnivariateDistribution
    design
    p
    function SampleSize{T<:Real}(design::AbstractBinaryTwoStageDesign, p::T)
        checkp(p)
        new(design, p)
    end
end


rand(d::SampleSize) = d.design |> interimsamplesize |> n1 -> Distributions.Binomial(n1, d.p) |> rand |> x1 -> samplesize(d.design, x1)


function pdf(d::SampleSize, n::Int)
    res = 0.0
    for x1 in 0:interimsamplesize(d.design)
        if samplesize(d.design, x1) == n
            res += Distributions.pdf(Distributions.Binomial(interimsamplesize(d.design), d.p), x1)
        end
    end
    return res
end


minimum(d::SampleSize) = interimsamplesize(d.design)

maximum(d::SampleSize) = maximum(samplesize(d.design))


function quantile(d::SampleSize, q::Real)
    nrange = minimum(d):maximum(d)
    res    = 0
    while Distributions.cdf(d, res + 1) <= q
        res += 1
    end
    return res
end


function mean(d::SampleSize)
    res = 0.0
    for n in minimum(d):maximum(d)
        res += pdf(d, n)*n
    end
    return(res)
end


function var(d::SampleSize)
    res = 0.0
    m   = mean(d)
    for n in minimum(d):maximum(d)
        res += pdf(d, n)*(n - m)^2
    end
    return(res)
end
