type SampleSize <: Distributions.DiscreteUnivariateDistribution
    design
    p
    function SampleSize{T<:Real}(design::AbstractBinaryTwoStageDesign, p::T)
        checkp(p)
        new(design, p)
    end
end


rand(d::SampleSize) = d.design |> getInterimSampleSize |> n1 -> Distributions.Binomial(n1, d.p) |> rand |> x1 -> getSampleSize(d.design, x1)


function pdf(d::SampleSize, n::Int)
    res = 0.0
    for x1 in 0:getInterimSampleSize(d.design)
        if getSampleSize(d.design, x1) == n
            res += Distributions.pdf(Distributions.Binomial(getInterimSampleSize(d.design), d.p), x1)
        end
    end
    return res
end


minimum(d::SampleSize) = getInterimSampleSize(d.design)

maximum(d::SampleSize) = maximum(getSampleSize(d.design))


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
