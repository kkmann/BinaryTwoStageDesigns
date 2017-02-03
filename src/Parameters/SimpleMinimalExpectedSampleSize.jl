type SimpleMinimalExpectedSampleSize{T_samplespace<:SampleSpace, T_reg<:Regularization, T_eff<:Efficacy} <: PointAlternative
    samplespace::T_samplespace
    p0
    p1
    alpha
    beta
    pess
    function SimpleMinimalExpectedSampleSize(
        samplespace,
        p0, p1,
        alpha, beta,
        pess
    )
        any(!([alpha; beta; p0; p1; pess] .>= 0.0)) ? throw(InexactError()) : nothing
        any(!([alpha; beta; p0; p1; pess] .<= 1.0)) ? throw(InexactError()) : nothing
        p0 >= p1 ? throw(InexactError()) : nothing
        new(samplespace, p0, p1, alpha, beta, pess)
    end
end
function SimpleMinimalExpectedSampleSize{T_samplespace<:SampleSpace}(
    samplespace::T_samplespace,
    p0, p1,
    alpha, beta,
    pess,
    T_reg, T_eff
)
    T_reg <: Regularization ? nothing : throw(InexactError()) # this must be possible in a more Julian way
    T_eff <: Efficacy ? nothing : throw(InexactError())
    SimpleMinimalExpectedSampleSize{T_samplespace, T_reg, T_eff}(samplespace, p0, p1, alpha, beta, pess)
end

maxsamplesize(params::SimpleMinimalExpectedSampleSize) = maxsamplesize(params.samplespace)
efficacy{T_samplespace, T_reg, T_eff}(params::SimpleMinimalExpectedSampleSize{T_samplespace, T_reg, T_eff}) = T_eff
regularization{T_samplespace, T_reg, T_eff}(params::SimpleMinimalExpectedSampleSize{T_samplespace, T_reg, T_eff}) = T_reg

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
    if efficacy(params) == StoppingForEfficacy
        cvals = [-Inf; 0:(nmax - 1); Inf]
    end
    if efficacy(params) == NoStoppingForEfficacy
        cvals = [0:(nmax - 1); Inf]
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
    # add optimality criterion
    @objective(m, Min,
        sum(dbinom(x1, n1, p1)*n*y[x1, n, c] for
            x1 in 0:n1,
            n  in nvals,
            c  in cvals
        )
    )
    return m, y
end
