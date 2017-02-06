# define types for options
abstract IsGroupSequential
type GroupSequential <: IsGroupSequential
end
type NotGroupSequential <: IsGroupSequential
end

abstract StoppingForEfficacy
type NoStoppingForEfficacy <: StoppingForEfficacy
end
type AllowStoppingForEfficacy <: StoppingForEfficacy
end

abstract HasMonotoneConditionalPower
type NoMonotoneConditionalPower <: HasMonotoneConditionalPower
end
type MonotoneConditionalPower <: HasMonotoneConditionalPower
end


type SimpleMinimalExpectedSampleSize{
        T_samplespace<:SampleSpace,
        T_gs<:IsGroupSequential,
        T_eff<:StoppingForEfficacy,
        T_mcp<:HasMonotoneConditionalPower
} <: PointAlternative
    samplespace::T_samplespace
    p0
    p1
    alpha
    beta
    pess
    mincondpower
    function SimpleMinimalExpectedSampleSize(
        samplespace,
        p0, p1,
        alpha, beta,
        pess,
        mincondpower
    )
        any(!([alpha; beta; p0; p1; pess; mincondpower] .>= 0.0)) ? throw(InexactError()) : nothing
        any(!([alpha; beta; p0; p1; pess; mincondpower] .<= 1.0)) ? throw(InexactError()) : nothing
        p0 >= p1 ? throw(InexactError()) : nothing
        new(samplespace, p0, p1, alpha, beta, pess, mincondpower)
    end
end
function SimpleMinimalExpectedSampleSize{T_samplespace<:SampleSpace}(
    samplespace::T_samplespace,
    p0, p1,
    alpha, beta,
    pess;
    minconditionalpower::Real = 0.0,
    GROUPSEQUENTIAL::Bool     = false,
    STOPPINGFOREFFICACY::Bool = true,
    MONOTONECONDITIONALPOWER::Bool = false
)
    T_gs = NotGroupSequential
    if GROUPSEQUENTIAL
        T_gs = GroupSequential
    end
    T_eff = AllowStoppingForEfficacy
    if !STOPPINGFOREFFICACY
        T_eff = NoStoppingForEfficacy
    end
    T_mcp = NoMonotoneConditionalPower
    if MONOTONECONDITIONALPOWER
        T_mcp = MonotoneConditionalPower
    end
    SimpleMinimalExpectedSampleSize{T_samplespace, T_gs, T_eff, T_mcp}(
        samplespace, p0, p1, alpha, beta, pess, minconditionalpower
    )
end

maxsamplesize(params::SimpleMinimalExpectedSampleSize) = maxsamplesize(params.samplespace)

allowsstoppingforefficacy{T_samplespace, T_gs, T_eff, T_mcp}(params::SimpleMinimalExpectedSampleSize{T_samplespace, T_gs, T_eff, T_mcp}) = T_eff == AllowStoppingForEfficacy ? true : false

isgroupsequential{T_samplespace, T_gs, T_eff, T_mcp}(params::SimpleMinimalExpectedSampleSize{T_samplespace, T_gs, T_eff, T_mcp}) = T_gs == GroupSequential ? true : false

hasmonotoneconditionalpower{T_samplespace, T_gs, T_eff, T_mcp}(params::SimpleMinimalExpectedSampleSize{T_samplespace, T_gs, T_eff, T_mcp}) = T_mcp == MonotoneConditionalPower ? true : false

minconditionalpower(params::SimpleMinimalExpectedSampleSize) = params.mincondpower

function score{T_P<:SimpleMinimalExpectedSampleSize}(design::AbstractBinaryTwoStageDesign, params::T_P)
    n = SampleSize(design, params.pess)
    return mean(n)
end

function _createProblem{T<:Integer}(
    n1::T,      # stage one sample size
    params::SimpleMinimalExpectedSampleSize
)
    ss = samplespace(params)
    nmax = maxsamplesize(ss)
    possible(n1, ss) ? nothing : throw(InexactError())
    a  = alpha(params)
    b  = beta(params)
    p0 = null(params)
    p1 = alternative(params)
    # define base problem
    m, y = _createBaseProblem(n1, params) # c.f. util.jl
    nvals = n1:nmax
    cvalsfinite = 0:(nmax - 1)
    if allowsstoppingforefficacy(params)
        cvals = [-Inf; cvalsfinite; Inf]
        cvalsinfinite = [-Inf; Inf]
    else
        cvals = [cvalsfinite; Inf]
        cvalsinfinite = [Inf]
    end

    # add type one error rate constraint
    @constraint(m,
        sum(dbinom(x1, n1, p0)*_cpr(x1, n1, n, c, p0)*y[x1, n, c] for
            x1 in 0:n1,
            n  in nvals,
            c  in cvals
        ) <= a
    )
    # add type two error rate constraint (power)
    @constraint(m,
        sum(dbinom(x1, n1, p1)*_cpr(x1, n1, n, c, p1)*y[x1, n, c] for
            x1 in 0:n1,
            n  in nvals,
            c  in cvals
        ) >= 1 - b
    )
    # add conditional type two error rate constraint (power)
    for x1 in 0:n1
        # ensure monotonicity if required
        if x1 > 0 & hasmonotoneconditionalpower(params)
            @constraint(m,
                sum(_cpr(x1, n1, n, c, p1)*(y[x1, n, c] - y[x1, n, c]) for
                    n  in nvals,
                    c  in cvals
                ) >= 0
            )
        end
        @constraint(m,
            sum(_cpr(x1, n1, n, c, p1)*y[x1, n, c] for
                n  in nvals,
                c  in cvalsfinite
            ) + sum(y[x1, n, c] for
                n  in nvals,
                c  in cvalsinfinite
            ) >= minconditionalpower(params)
        )
    end
    # add optimality criterion
    @objective(m, Min,
        sum(dbinom(x1, n1, params.pess)*n*y[x1, n, c] for
            x1 in 0:n1,
            n  in nvals,
            c  in cvals
        )
    )
    return m, y
end

# TODO: expand and make compulsory!
function _isfeasible(design::BinaryTwoStageDesign, params::SimpleMinimalExpectedSampleSize)
    all(power.(design, linspace(0, null(params))) .<= alpha(params) + .001) ? nothing : throw(InexactError())
    all(power.(design, linspace(alternative(params), 1)) .>= 1 - beta(params) - .001) ? nothing : throw(InexactError())
    return true
end
