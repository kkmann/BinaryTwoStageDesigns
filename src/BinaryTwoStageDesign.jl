"""
Basic immutable class for two-stage designs with binary outcome variables.
"""
immutable BinaryTwoStageDesign  # need T2 must be able to hold +/- infinity
    n
    c
    function BinaryTwoStageDesign{T1<:Integer, T2<:Real}(n::Vector{T1}, c::Vector{T2})
        @assert all(n .>= length(n) - 1)
        @assert length(n) == length(c)
        new(n, c)
    end
end

function show(io::IO, object::BinaryTwoStageDesign)
    println("Binary two-stage design:\n=================================")
    println("      x1  n(x1) c(x1)")
    println("    ----- ----- -----")
    for x1 in 0:interimSampleSize(object)
        println(@sprintf "    %5i %5i %5i" x1 finalSampleSize(object, x1) rejectionBoundary(object, x1))
    end
    println("    ----- ----- -----")
end

function convert(::DataFrames.DataFrame, design::BinaryTwoStageDesign)
    df = DataFrames.DataFrame(
        x1 = 0:interimSampleSize(design),
        n  = finalSampleSize(design),
        c  = rejectionBoundary(design)
    )
    return df
end

function interimSampleSize(design::BinaryTwoStageDesign)
    return length(design.n) - 1
end

function finalSampleSize(design::BinaryTwoStageDesign)
    return design.n
end
function finalSampleSize{T1<:Integer}(design::BinaryTwoStageDesign, x1::T1)
    @assert 0 <= x1
    @assert x1 <= interimSampleSize(design)
    return design.n[x1 + 1] # array indexing starts at 1, not 0
end

function rejectionBoundary(design::BinaryTwoStageDesign)
    return design.c
end
function rejectionBoundary{T1<:Integer}(design::BinaryTwoStageDesign, x1::T1)
    @assert 0 <= x1
    @assert x1 <= interimSampleSize(design)
    return design.c[x1 + 1] # array indexing starts at 1, not 0
end



"""
Compute the conditional rejection probability given number of stage one
responses (x1) and response probability p.
"""
function conditionalProbabilityToReject{T1<:Integer, T2<:Real}(design::BinaryTwoStageDesign, x1::Vector{T1}, p::T2)
    @assert 0.0 <= p
    @assert p   <= 1.0
    n1 = interimSampleSize(design)
    @assert 0   <= x1
    @assert x1  <= n1
    c  = rejectionBoundary(design, x1)
    n  = finalSampleSize(design, x1)
    return _cpr(x1, n1, n, c, p) # defined in util.jl
end

"""
Compute the rejection probability given response probability p.
"""
function probabilityToReject{T<:Real}(design::BinaryTwoStageDesign, p::T)
    n1      = interimSampleSize(design)
    X1      = Distributions.Binomial(n1, p) # stage one responses
    x1range = collect(0:n1)
    return vecdot(Distributions.pdf(X1, x1range), conditionalProbabilityToReject.(design, x1, p))
end

"""
Create random variable for the final sample size
"""
type FinalSampleSize <: Distributions.DiscreteUnivariateDistribution
    design
    p
    function FinalSampleSize{T<:Real}(design::BinaryTwoStageDesign, p::T)
        @assert 0.0 <= p
        @assert p <= 1.0
        new(design, p)
    end
end

function rand(d::FinalSampleSize)
    return d.design |> interimSampleSize |> n1 -> Distributions.Binomial(n1, p) |> rand |> x1 -> finalSampleSize(d.design, x1)
end

function pdf(d::FinalSampleSize, n::Real)
    res = 0.0
    for x1 in 0:interimSampleSize(d.design)
        if finalSampleSize(d.design, x1) == n
            res += Distributions.pdf(Binomial(interimSampleSize(d.design), p), x1)
        end
    end
    return res
end

function minimum(d::FinalSampleSize)
    return interimSampleSize(d.design)
end

function maximum(d::FinalSampleSize)
    return maximum(finalSampleSize(d.design))
end

function cdf(d::FinalSampleSize, n::Real)
    nrange = minimum(d):n
    res    = 0.0
    for n_ in nrange
        res += Distributions.pdf(d, n_)
    end
    return res
end

function quantile(d::FinalSampleSize, q::Real)
    nrange = minimum(d):maximum(d)
    res    = 0
    while Distributions.cdf(d, res + 1) <= q
        res += 1
    end
    return res
end



"""
Check inputs (n, x) for reachability with given design

This method returns "true" iff the combination (n, x) is reachable with the
given design.
"""
function reachable{T<:Integer}(design::BinaryTwoStageDesign, n::T, x::T)
    if x < 0 | x > n | !(n in finalSampleSize(design))
        return false
    else
        res = false
        n1 = interimSampleSize(design)
        if n == n1 # only first stage
            if finalSampleSize(design, x) == n1
                res = true
            end
        else # also second stage
            for x1 in 0:(min(x, n1))
                if (finalSampleSize(design, x1) == n) & (x - x1 <= finalSampleSize(design, x1) - n1)
                    res = true
                end
            end
        end
    end
    return res
end

"""
Check inputs (x1, x2) for reachability with given design

This method returns "true" iff the combination (x1, x2) is reachable with the
given design.
"""
function reachableResponses{T<:Integer}(design::BinaryTwoStageDesign, x1::T, x2::T)
    n1 = interimSampleSize(design)
    if (x1 < 0) | (x1 > n1)
        return false
    end
    if (x2 > finalSampleSize(design, x1) - n1) | x2 < 0
        return false
    else
        return true
    end
end

"""
Probability mass function of a given design
"""
function pmf_x1x2(design::BinaryTwoStageDesign, x1::Int64, x2::Int64, p::Float64)
    n1::Int64 = length(design.n) - 1
    n::Int64 = design.n[x1 + 1]
    if isPossible_x1x2(design, x1, x2)[1]
        res = p^(x1 + x2)*(1 - p)^(n - x1 - x2)*binomial(BigInt(n1), BigInt(x1))*binomial(BigInt(n - n1), BigInt(x2))
    else
        # impossible has PMF 0.0
        res = 0.0
    end
    return(res)
end

"""
Test

Obtain a test decision for a given design and observed stage one and stage two outcomes. The
method returns 'true' whenever the null hypothesis is rejected and 'false' otherwise. Inputs are
checked for consistency.
"""
function test(design::BinaryTwoStageDesign, x1::Int64, x2::Int64 = 0)
    n1::Int64 = length(design.n) - 1
    n::Int64 = design.n[x1 + 1]
    c::Float64 = design.c[x1 + 1]
    if (x1 < 0) | (x1 > n1) > ~isfinite(x1)
        return(NaN)
    end
    if (x2 < 0) | (x2 > n - n1) | ~isfinite(x2)
        return(NaN)
    end
    if x1 + x2 > c
        return true
    else
        return false
    end
end
function test(design::BinaryTwoStageDesign, x1::Vector{Int64}, x2::Vector{Int64})
    # vectorized version
    m = length(x1)
    res = zeros(Float64, m)
    for i in 1:m
        res[i] = test(design, x1[i], x2[i])
    end
    return(res)
end
