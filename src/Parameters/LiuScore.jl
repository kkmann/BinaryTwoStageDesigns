type LiuScore{
        T_samplespace<:SampleSpace,
        T_gs<:IsGroupSequential,
        T_eff<:StoppingForEfficacy,
        T_mcp<:HasMonotoneConditionalPower
} <: PointAlternative
    samplespace::T_samplespace
    p0
    alpha
    beta
    p1
    mincondpower
    fs
    fp
    function LiuScore(
        samplespace,
        p0, p1,
        alpha, beta,
        mincondpower,
        fs, fp
    )
        any(!([alpha; beta; p0; p1; mincondpower] .>= 0.0)) ? throw(InexactError()) : nothing
        any(!([alpha; beta; p0; p1; mincondpower] .<= 1.0)) ? throw(InexactError()) : nothing
        (fs >= 0) & (fs < 1) ? nothing : throw(InexactError())
        (fp >= 0) & (fp < 1) ? nothing : throw(InexactError())
        new(samplespace, p0, alpha, beta, p1, mincondpower, fs, fp)
    end
end
function LiuScore{T_samplespace<:SampleSpace}(
    samplespace::T_samplespace,
    p0, p1, alpha, beta,
    fs, fp;
    minconditionalpower::Real      = 0.0,
    GROUPSEQUENTIAL::Bool          = false,
    STOPPINGFOREFFICACY::Bool      = true,
    MONOTONECONDITIONALPOWER::Bool = true
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
    LiuScore{T_samplespace, T_gs, T_eff, T_mcp}(
        samplespace, p0, p1, alpha, beta, minconditionalpower, fs, fp
    )
end

maxsamplesize(params::LiuScore) = maxsamplesize(params.samplespace)

allowsstoppingforefficacy{T_samplespace, T_gs, T_eff, T_mcp}(params::LiuScore{T_samplespace, T_gs, T_eff, T_mcp}) = T_eff == AllowStoppingForEfficacy ? true : false

isgroupsequential{T_samplespace, T_gs, T_eff, T_mcp}(params::LiuScore{T_samplespace, T_gs, T_eff, T_mcp}) = T_gs == GroupSequential ? true : false

hasmonotoneconditionalpower{T_samplespace, T_gs, T_eff, T_mcp}(params::LiuScore{T_samplespace, T_gs, T_eff, T_mcp}) = T_mcp == MonotoneConditionalPower ? true : false

minconditionalpower(params::LiuScore) = params.mincondpower

beta(params::LiuScore) = params.beta

function score{T_P<:LiuScore}(design::AbstractBinaryTwoStageDesign, params::T_P)
    # ROS
    powerreq = 1 - beta(params)
    n1 = interimsamplesize(design)
    nreq__(p, powerreq) = nreq_(p, powerreq, params.p0, params.alpha)
    ros = 0.0
    for x1 in 0:n1
        ros += dbinom(x1, n1, params.p1)*max(0.0, samplesize(design, x1)/nreq__(params.p1, powerreq) - 1)
    end
    ros = 100/(1 - params.fs)*ros
    # rup
    rup = max(0, nreq__(params.p1, powerreq) - nreq__(params.p1, power(design, params.p1))) / (nreq__(params.p1, powerreq) - nreq__(p, (1 - fp)*powerreq))
    return ros + rup
end


function _createProblem{T<:Integer}(
    n1::T,      # stage one sample size
    params::LiuScore
)
    ss = samplespace(params)
    nmax = maxsamplesize(ss)
    possible(n1, ss) ? nothing : throw(InexactError())
    a  = alpha(params)
    b  = beta(params)
    p0 = null(params)
    p1 = params.p1
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
    nreq__(p, powerreq) = nreq_(p, powerreq, p0, a)
    # construct expressions for ROS, upper bound necessary to guarantee finteness!
    @expression(m, ros,
        100/(1 - params.fs)*sum(dbinom(x1, n1, p1)*min(999999, max(0, n / nreq__(p1, 1 - b) - 1 ))*y[x1, n, c] for
            x1 in 0:n1,
            n in nvals,
            c in cvals
        )
    )
    # construct expressions for RUP
    # piecewise linear approximation of h_p(power)
    # BULLSHIT USE SOS2 set to model arbitrary piecewise linear!
    @expression(m, power,
        sum(dbinom(x1, n1, p1)*_cpr(x1, n1, n, c, p1)*y[x1, n, c] for
            x1 in 0:n1,
            n  in nvals,
            c  in cvals
        )
    )
    pivots = collect(linspace(0, 1, 11))
    @variable(m, 0 <= lambda[piv = pivots] <= 1)
    addSOS2(m, [lambda[piv] for piv in pivots])
    @constraint(m, sum(lambda[piv] for piv in pivots) == 1)
    @constraint(m, sum(lambda[piv]*piv for piv in pivots) == power) # defines lambdas!
    f(power) = 100*max(0, nreq__(p1, 1 - b) - nreq__(p1, power)) / (nreq__(p1, 1 - b) - nreq__(p1, (1 - params.fp)*(1 - b)))
    @variable(m, rup)
    @constraint(m, sum(lambda[piv]*f(piv) for piv in pivots) == rup)
    @objective(m, Min,
        rup + ros
    )
    return m, y
end

function _isfeasible(design::BinaryTwoStageDesign, params::LiuScore)
    return true
end


# utility
function nreq_(p, power, p0, alpha)
    if power > alpha
        return (1 - p)*p*( (qnorm(1 - alpha) + qnorm(power)) / (p0 - p) )^2
    else
        return 0.0
    end
end
