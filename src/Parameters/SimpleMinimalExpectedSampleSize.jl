type SimpleMinimalExpectedSampleSize{T_samplespace<:SampleSpace} <: PointAlternative
    samplespace::T_samplespace
    p0
    p1
    alpha
    beta
    pess
    minconditionalpower
    MONOTONECONDITIONALPOWER
    minstoppingforfutility
    function SimpleMinimalExpectedSampleSize(
        samplespace,
        p0, p1,
        alpha, beta,
        pess,
        minconditionalpower,
        MONOTONECONDITIONALPOWER,
        minstoppingforfutility,
    )
        @checkprob p0 p1 alpha beta pess minconditionalpower minstoppingforfutility
        p0 >= p1 ? error("p0 must be smaller than p1") : nothing
        new(samplespace, p0, p1, alpha, beta, pess, minconditionalpower, MONOTONECONDITIONALPOWER, minstoppingforfutility)
    end
end
function SimpleMinimalExpectedSampleSize{T_samplespace<:SampleSpace}(
    samplespace::T_samplespace,
    p0, p1,
    alpha, beta,
    pess;
    minstoppingforfutility::Real   = 0.0,
    minconditionalpower::Real      = 0.0,
    MONOTONECONDITIONALPOWER::Bool = true
)
    SimpleMinimalExpectedSampleSize{T_samplespace}(
        samplespace, p0, p1, alpha, beta, pess, minconditionalpower, MONOTONECONDITIONALPOWER, minstoppingforfutility
    )
end

maxsamplesize(params::SimpleMinimalExpectedSampleSize) = maxsamplesize(params.samplespace)
isgroupsequential{T_samplespace}(params::SimpleMinimalExpectedSampleSize{T_samplespace}) = isgroupsequential(params.ss)
hasmonotoneconditionalpower{T_samplespace}(params::SimpleMinimalExpectedSampleSize{T_samplespace}) = params.MONOTONECONDITIONALPOWER
minconditionalpower(params::SimpleMinimalExpectedSampleSize) = params.minconditionalpower


function score{T_P<:SimpleMinimalExpectedSampleSize}(design::BinaryTwoStageDesign, params::T_P)
    n = SampleSize(design, params.pess)
    return mean(n)
end

function completemodel{T<:Integer}(ipm::IPModel, params::SimpleMinimalExpectedSampleSize, n1::T)
    ss = samplespace(params)
    !possible(n1, ss) ? warn("completemodel(EB): n1 and sample space incompatible") : nothing
    # extract ip model
    m             = ipm.m
    y             = ipm.y
    nvals         = ipm.nvals
    cvals         = ipm.cvals
    cvalsfinite   = ipm.cvalsfinite
    cvalsinfinite = [-Inf; Inf]
    # extract other parameters
    nmax          = maxsamplesize(ss, n1)
    p0            = null(params)
    p1            = params.p1
    # add type one error rate constraint
    JuMP.@constraint(m,
        sum(dbinom(x1, n1, p0)*_cpr(x1, n1, n, c, p0)*y[x1, n, c] for
            x1 in 0:n1,
            n  in nvals,
            c  in cvals
        ) <= alpha(params)
    )
    # add type two error rate constraint (power)
    JuMP.@constraint(m,
        sum(dbinom(x1, n1, p1)*_cpr(x1, n1, n, c, p1)*y[x1, n, c] for
            x1 in 0:n1,
            n  in nvals,
            c  in cvals
        ) >= 1 - params.beta
    )
    # add conditional type two error rate constraint (power)
    for x1 in 0:n1
        # ensure monotonicity if required
        if x1 >= 1 & hasmonotoneconditionalpower(params)
            JuMP.@constraint(m,
                sum(_cpr(x1, n1, n, c, p1)*y[x1, n, c] - _cpr(x1 - 1, n1, n, c, p1)*y[x1 - 1, n, c] for
                    n  in nvals, c in cvals
                ) >= 0
            )
        end
        JuMP.@constraint(m,
            sum(_cpr(x1, n1, n, c, p1)*y[x1, n, c] for
                n  in nvals,
                c  in cvalsfinite
            ) + sum(y[x1, n, c] for
                n  in nvals,
                c  in cvalsinfinite
            ) >= minconditionalpower(params)
        )
    end
    # add constraint for minimal stopping-for-futility probability
    JuMP.@constraint(m,
        sum(dbinom(x1, n1, p0)*y[x1, n1, Inf] for
            x1 in 0:n1
        ) >= params.minstoppingforfutility
    )
    # add optimality criterion
    JuMP.@objective(m, Min,
        sum(dbinom(x1, n1, params.pess)*n*y[x1, n, c] for
            x1 in 0:n1,
            n  in nvals,
            c  in cvals
        )
    )
    return true
end

function _isfeasible(design::BinaryTwoStageDesign, params::SimpleMinimalExpectedSampleSize)
    all(power.(design, linspace(0, null(params))) .<= alpha(params) + .001) ? nothing : throw(InexactError())
    all(power.(design, linspace(alternative(params), 1)) .>= 1 - beta(params) - .001) ? nothing : throw(InexactError())
    return true
end
