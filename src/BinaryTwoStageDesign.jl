abstract AbstractBinaryTwoStageDesign


# these two enable array broadcasting of all methods!
Base.size(::AbstractBinaryTwoStageDesign) = ()
Base.length(design::AbstractBinaryTwoStageDesign) = 1
Base.start(design::AbstractBinaryTwoStageDesign) = true
Base.done(design::AbstractBinaryTwoStageDesign, state) = true
Base.getindex(design::AbstractBinaryTwoStageDesign, i) = design


immutable BinaryTwoStageDesign{T1<:Integer, T2<:Real, PType<:Parameters} <: AbstractBinaryTwoStageDesign
    n::Vector{T1}
    c::Vector{T2} # rejected if x1 + x2 > c(x1), must allow real to hold Infinity
    params::PType
    function BinaryTwoStageDesign(n, c)
        any(n .< (length(n) - 1)) ? throw(InexactError()) : nothing
        length(n) != length(c) ? throw(InexactError()) : nothing
        new(n, c, NoParameters())
    end
    function BinaryTwoStageDesign(n, c, params)
        any(n .< (length(n) - 1)) ? throw(InexactError()) : nothing
        length(n) != length(c) ? throw(InexactError()) : nothing
        new(n, c, params)
    end
end # BinaryTwoStageDesign
BinaryTwoStageDesign{T1<:Integer, T2<:Real}(
    n::Vector{T1},
    c::Vector{T2}
) = BinaryTwoStageDesign{T1, T2, NoParameters}(n, c)
BinaryTwoStageDesign{T1<:Integer, T2<:Real, PType<:Parameters}(
    n::Vector{T1},
    c::Vector{T2},
    params::PType
) = BinaryTwoStageDesign{T1, T2, PType}(n, c, params)

convert(::Type{DataFrames.DataFrame}, design::AbstractBinaryTwoStageDesign) =
    DataFrames.DataFrame(
        x1 = 0:interimsamplesize(design),
        n  = samplesize(design),
        c  = criticalvalue(design)
    )


interimsamplesize(design::AbstractBinaryTwoStageDesign) = length(design.n) - 1
parameters(design::AbstractBinaryTwoStageDesign) = design.params

samplesize(design::AbstractBinaryTwoStageDesign) = design.n
samplesize{T<:Integer}(design::AbstractBinaryTwoStageDesign, x1::T) =
    (checkx1(x1, design); design.n[x1 + 1])


criticalvalue(design::AbstractBinaryTwoStageDesign) = design.c
criticalvalue{T<:Integer}(design::AbstractBinaryTwoStageDesign, x1::T) =
    (checkx1(x1, design); design.c[x1 + 1])


function pdf{T1<:Integer, T2<:Real}(
    design::AbstractBinaryTwoStageDesign,
    x1::T1, x2::T1, p::T2
)
    checkp(p)
    try
        checkx1x2(x1, x2, design)
    catch
        return 0.0
    end
    n1 = interimsamplesize(design)
    n  = samplesize(design, x1)
    try
        return p^(x1 + x2)*(1 - p)^(n - x1 - x2)*binomial(n1, x1)*binomial(n - n1, x2) # is 0 if impossible
    catch
        return p^(x1 + x2)*(1 - p)^(n - x1 - x2)*binomial(BigInt(n1), BigInt(x1))*binomial(BigInt(n - n1), BigInt(x2)) # is 0 if impossible
    end

end


function power{T1<:Integer, T2<:Real}(
    design::AbstractBinaryTwoStageDesign, x1::T1, p::T2
)
    checkp(p)
    checkx1(x1, design)
    res = _cpr(
        x1,
        interimsamplesize(design),
        samplesize(design, x1),
        criticalvalue(design, x1),
        p
    )
    return min(1, max(0, res)) # guarantee bounds!
end
function power{T<:Real}(design::AbstractBinaryTwoStageDesign, p::T)
    checkp(p)
    n1      = interimsamplesize(design)
    X1      = Distributions.Binomial(n1, p) # stage one responses
    x1range = collect(0:n1)
    return min(1, max(0, vecdot(Distributions.pdf(X1, x1range), power.(design, x1range, p))))
end

function expectedpower(design::AbstractBinaryTwoStageDesign, params::VagueAlternative, x1)
    pmcrv    = mcrv(params)
    phi(p)   = prior(params, p)
    n1       = interimsamplesize(design)
    z        = quadgk(p -> phi(p) * Distributions.pdf(Distributions.Binomial(n1, p), x1), pmcrv, 1)[1]
    omega(p) = p < pmcrv ? 0 : phi(p) * Distributions.pdf(Distributions.Binomial(n1, p), x1) / z # conditional prior for stage 2
    return quadgk(p -> power(design, x1, p) * omega(p), pmcrv, 1)[1]
end
function expectedpower(design::AbstractBinaryTwoStageDesign, params::VagueAlternative)
    pmcrv  = mcrv(params)
    phi(p) = prior(params, p)
    n1     = interimsamplesize(design)
    z      = quadgk(phi, pmcrv, 1)[1]
    omega(p) = p < pmcrv ? 0 : phi(p) / z # conditional prior given effect
    return quadgk(p -> power(design, p) * omega(p), pmcrv, 1)[1]
end


function stoppingforfutility{T<:Real}(design::AbstractBinaryTwoStageDesign, p::T)
    checkp(p)
    n1  = interimsamplesize(design)
    X1  = Distributions.Binomial(n1, p) # stage one responses
    res = 0.0
    c   = criticalvalue(design)
    for x1 in 0:n1
        if c[x1 + 1] == Inf
            res += Distributions.pdf(X1, x1)
        end
    end
    return res
end


function test{T<:Integer}(design::AbstractBinaryTwoStageDesign, x1::T, x2::T)::Bool
    checkx1x2(x1, x2, design)
    return x1 + x2 > criticalvalue(design, x1)
end


function simulate{T1<:Integer, T2<:Real}(design::AbstractBinaryTwoStageDesign, p::T2, nsim::T1)
    x2    = SharedArray(Int, nsim)
    n     = SharedArray(Int, nsim)
    c     = SharedArray(Float64, nsim)
    rej   = SharedArray(Bool, nsim)
    n1    = interimsamplesize(design)
    rv_x1 = Distributions.Binomial(n1, p)
    x1    = SharedArray(Int, nsim)
    @sync @parallel for i in 1:nsim
        x1[i]  = rand(rv_x1)
        n2     = samplesize(design, x1[i]) - n1
        n[i]   = n2 + n1
        c[i]   = criticalvalue(design, x1[i])
        x2[i]  = rand(Distributions.Binomial(n2, p))
        rej[i] = test(design, x1[i], x2[i])
    end
    return DataFrames.DataFrame(
        x1 = convert(Vector{Int}, x1),
        n  = convert(Vector{Int}, n),
        c  = convert(Vector{Float64}, c),
        x2 = convert(Vector{Int}, x2),
        rejectedH0 = convert(Vector{Bool}, rej)
    )
end


score(design::AbstractBinaryTwoStageDesign, params::Parameters)::Real = error("not implemented")
score(design::AbstractBinaryTwoStageDesign)::Real = score(design, parameters(design))


function jeffreysprior(design::AbstractBinaryTwoStageDesign)
    function sqrtfi(p::Float64)::Float64
        res  = 0.0
        n1   = interimsamplesize(design)
        nmax = maximum(samplesize(design))
        supp = support(design)
        for i in 1:size(supp, 1)
            x1, x2 = supp[i, 1], supp[i, 2]
            x   = x1 + x2
            n   = samplesize(design, x1)
            res += binomial(BigInt(n1), BigInt(x1))*binomial(BigInt(n - n1), BigInt(x2))*BigFloat(p^x*(1 - p)^(n - x)*(x/p - (n - x)/(1 - p))^2)
        end
        return sqrt(res)
    end
    z = quadgk(sqrtfi, 0, 1, abstol = 0.001)[1] # exact integration from 0 to 1 is expensive!
    function prior{T<:Real}(p::T)::Float64
        checkp(p)
        return sqrtfi(p)/z
    end
    return prior
end



# utility functions
function checkx1{T<:Integer}(x1::T, design::AbstractBinaryTwoStageDesign)
    x1 < 0 ? throw(InexactError("x1 must be non-negative")) : nothing
    x1 > interimsamplesize(design) ? throw(InexactError("x1 smaller or equal to n1")) : nothing
end

function checkx1x2{T<:Integer}(x1::T, x2::T, design::AbstractBinaryTwoStageDesign)
    checkx1(x1, design)
    x2 < 0 ? throw(InexactError("x2 must be non-negative")) : nothing
    n1 = interimsamplesize(design)
    n2 = samplesize(design, x1) - n1
    x2 > n2 ? throw(InexactError("x2 must be smaller or equal to n2")) : nothing
end

function ispossible{T<:Integer}(design::AbstractBinaryTwoStageDesign, x1::T, x2::T)
    res = x1 < 0 ? false : true
    res = x1 > interimsamplesize(design) ? false : true
    res = x2 < 0 ? false : true
    res = x2 > samplesize(design, x1) - interimsamplesize(design) ? false : true
end

function support(design::AbstractBinaryTwoStageDesign)
    n1     = interimsamplesize(design)
    nmax   = maximum(samplesize(design))
    return [[x1, x2] for x1 in 0:n1, x2 in 0:(nmax - n1) if
        (x2 <= samplesize(design, x1) - n1) & ((samplesize(design, x1) > n1) | (x2 == 0))
    ] |> x-> hcat(x...)'
end

checkp{T<:Real}(p::T) = (0.0 > p) & (p > 1.0) ? throw(InexactError("p must be in [0, 1]")) : nothing

function _cpr(x1, n1, n, c, p)
    # conditional probability to reject
    if x1 > c
        return(1.0)
    end
    if n - n1 + x1 <= c
        return(0.0)
    end # TODO: simplify
    return 1 - Distributions.cdf(Distributions.Binomial(n - n1, p), convert(Int64, c - x1))
end

function _cprprior(x1, n1, n, c, p0, prior)
    if x1 > c
        return 1.0
    end
    if n - n1 + x1 <= c
        return 0.0
    end
    z = QuadGK.quadgk(p -> prior(p)*dbinom(c - x1, n - n1, p), p0, 1.0, abstol = 1e-4)[1]
    cposterior(p) = prior(p)*dbinom(c - x1, n - n1, p)/z
    return QuadGK.quadgk(p -> cposterior(p)*_cpr(x1, n1, n, c, p), p0, 1.0, abstol = 1e-4)[1]
end
