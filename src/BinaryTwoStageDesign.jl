abstract AbstractBinaryTwoStageDesign

# these two enable array broadcasting of all methods!
Base.size(::AbstractBinaryTwoStageDesign) = ()
Base.getindex(design::AbstractBinaryTwoStageDesign, i) = design

immutable BinaryTwoStageDesign <: AbstractBinaryTwoStageDesign # need T2 must be able to hold +/- infinity
    n
    c # rejected if x1 + x2 > c(x1)
    function BinaryTwoStageDesign{T1<:Integer, T2<:Real}(n::Vector{T1}, c::Vector{T2})
        # @assert all(n .>= length(n) - 1)
        if any(n .< length(n) - 1)
            println(DataFrames.DataFrame(n = n, c = c))
            error("n must be greater or equal to n1")
        end
        @assert length(n) == length(c)
        new(n, c)
    end
end

function show(io::IO, object::AbstractBinaryTwoStageDesign)
    println("Binary two-stage design:\n=================================")
    println("      x1  n(x1) c(x1)")
    println("    ----- ----- -----")
    for x1 in 0:getInterimSampleSize(object)
        println(@sprintf "    %5i %5i %5i" x1 getSampleSize(object, x1) getRejectionBoundary(object, x1))
    end
    println("    ----- ----- -----")
end

function convert(::Type{DataFrames.DataFrame}, design::AbstractBinaryTwoStageDesign)
    df = DataFrames.DataFrame(
        x1 = 0:getInterimSampleSize(design),
        n  = getSampleSize(design),
        c  = getRejectionBoundary(design)
    )
    return df
end

function getInterimSampleSize(design::AbstractBinaryTwoStageDesign)
    return length(design.n) - 1
end

function getSampleSize(design::AbstractBinaryTwoStageDesign)
    return design.n
end
function getSampleSize{T1<:Integer}(design::AbstractBinaryTwoStageDesign, x1::T1)
    @assert 0 <= x1
    @assert x1 <= getInterimSampleSize(design)
    return design.n[x1 + 1] # array indexing starts at 1, not 0
end

function getRejectionBoundary(design::AbstractBinaryTwoStageDesign) # TODO: rename: critical value
    return design.c
end
function getRejectionBoundary{T1<:Integer}(design::AbstractBinaryTwoStageDesign, x1::T1)
    @assert 0 <= x1
    @assert x1 <= getInterimSampleSize(design)
    return design.c[x1 + 1] # array indexing starts at 1, not 0
end

function probability{T1<:Integer, T2<:Real}(design::AbstractBinaryTwoStageDesign, x1::T1, x2::T1, p::T2)::Float64
    if (x1 < 0) | (x2 < 0)
        return 0.0
    end
    n1 = getInterimSampleSize(design)
    n  = getSampleSize(design, x1)
    try
        return p^(x1 + x2)*(1 - p)^(n - x1 - x2)*binomial(n1, x1)*binomial(n - n1, x2) # is 0 if impossible
    catch
        return p^(x1 + x2)*(1 - p)^(n - x1 - x2)*binomial(BigInt(n1), BigInt(x1))*binomial(BigInt(n - n1), BigInt(x2)) # is 0 if impossible
    end

end

"""
Compute the conditional rejection probability given number of stage one
responses (x1) and response probability p.
"""
function conditionalProbabilityToReject{T1<:Integer, T2<:Real}(design::AbstractBinaryTwoStageDesign, x1::T1, p::T2)
    @assert 0.0 <= p
    @assert p   <= 1.0
    n1 = getInterimSampleSize(design)
    @assert 0   <= x1
    @assert x1  <= n1
    c  = getRejectionBoundary(design, x1)
    n  = getSampleSize(design, x1)
    return _cpr(x1, n1, n, c, p) # defined in util.jl
end

"""
Compute the rejection probability given response probability p.
"""
function probabilityToReject{T<:Real}(design::AbstractBinaryTwoStageDesign, p::T)
    n1      = getInterimSampleSize(design)
    X1      = Distributions.Binomial(n1, p) # stage one responses
    x1range = collect(0:n1)
    return vecdot(Distributions.pdf(X1, x1range), conditionalProbabilityToReject.(design, x1range, p))
end

function test{T<:Integer}(design::AbstractBinaryTwoStageDesign, x1::T, x2::T)::Bool
    supp = _support(design)
    if !_isInSupport(supp, x1, x2)
        error("observations are not in designs support")
    else
        return x1 + x2 > getRejectionBoundary(design, x1) ? true : false
    end
end

function simulate{T1<:Integer, T2<:Real}(design::AbstractBinaryTwoStageDesign, p::T2, nsim::T1)
    x2    = SharedArray(Int, nsim)
    n     = SharedArray(Int, nsim)
    c     = SharedArray(Float64, nsim)
    rej   = SharedArray(Bool, nsim)
    n1    = getInterimSampleSize(design)
    rv_x1 = Distributions.Binomial(n1, p)
    x1    = SharedArray(Int, nsim)
    @sync @parallel for i in 1:nsim
        x1[i]  = rand(rv_x1)
        n2     = getSampleSize(design, x1[i]) - n1
        n[i]   = n2 + n1
        c[i]   = getRejectionBoundary(design, x1[i])
        x2[i]  = rand(Distributions.Binomial(n2, p))
        rej[i] = test(design, x1[i], x2[i])
    end
    # return x1, n, c, x2, rej
    return DataFrames.DataFrame(
        x1 = convert(Vector{Int}, x1),
        n  = convert(Vector{Int}, n),
        c  = convert(Vector{Float64}, c),
        x2 = convert(Vector{Int}, x2),
        rejectedH0 = convert(Vector{Bool}, rej)
    )
end

"""
Create random variable for the final sample size
"""
type SampleSize <: Distributions.DiscreteUnivariateDistribution
    design
    p
    function SampleSize{T<:Real}(design::AbstractBinaryTwoStageDesign, p::T)
        @assert 0.0 <= p
        @assert p <= 1.0
        new(design, p)
    end
end

function rand(d::SampleSize)
    return d.design |> getInterimSampleSize |> n1 -> Distributions.Binomial(n1, d.p) |> rand |> x1 -> getSampleSize(d.design, x1)
end

function pdf(d::SampleSize, n::Int)
    res = 0.0
    for x1 in 0:getInterimSampleSize(d.design)
        if getSampleSize(d.design, x1) == n
            res += Distributions.pdf(Distributions.Binomial(getInterimSampleSize(d.design), d.p), x1)
        end
    end
    return res
end

function minimum(d::SampleSize)
    return getInterimSampleSize(d.design)
end

function maximum(d::SampleSize)
    return maximum(getSampleSize(d.design))
end

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
