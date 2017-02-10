type MinimalMinimalExpectedSampleSize{T_samplespace<:SampleSpace} <: PointAlternative
    samplespace::T_samplespace
    p0
    p1
    alpha
    beta
    pess
    function MinimalMinimalExpectedSampleSize(
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
function MinimalMinimalExpectedSampleSize{T_samplespace<:SampleSpace}(
    samplespace::T_samplespace,
    p0, p1,
    alpha, beta,
    pess
)
    MinimalMinimalExpectedSampleSize{T_samplespace}(
        samplespace, p0, p1, alpha, beta, pess
    )
end

maxsamplesize(params::MinimalMinimalExpectedSampleSize) = maxsamplesize(params.samplespace)

allowsstoppingforefficacy{T_samplespace}(params::MinimalMinimalExpectedSampleSize{T_samplespace}) = true

isgroupsequential{T_samplespace}(params::MinimalMinimalExpectedSampleSize{T_samplespace}) = false

function score{T_P<:MinimalMinimalExpectedSampleSize}(design::AbstractBinaryTwoStageDesign, params::T_P)
    n = SampleSize(design, params.pess)
    return mean(n)
end

function _createProblem{T<:Integer}(
    n1::T,      # stage one sample size
    params::MinimalMinimalExpectedSampleSize
)
    ss   = samplespace(params)
    nmax = maxsamplesize(ss)
    possible(n1, ss) ? nothing : throw(InexactError())
    a  = alpha(params)
    b  = beta(params)
    p0 = null(params)
    p1 = alternative(params)
    # define base problem
    !possible(n1, ss) ? throw(InexactError()) : nothing
    nmax = maxsamplesize(ss)
    nvals = n1:nmax
    cvalsfinite = 0:(nmax - 1)
    cvals = [-Inf; cvalsfinite; Inf]
    cvalsinfinite = [-Inf; Inf]
    m = Model()
    @variable(m, # indicator variables y[x1, n, c] == 1 iff n(x1) = n, c(x1) = c
        y[x1 = 0:n1, n = nvals, c = cvals],
        Bin
    )
    for x1 in 0:n1
        @constraint(m, # functional constraint: exactly one non-zero entry in (y) per x1
            sum(y[x1, n, c] for n in nvals, c in cvals) == 1
        )
        for n in nvals
            for c in cvals
                if isfinite(c)
                    if c >= n
                        @constraint(m, y[x1, n, c] == 0)
                    end
                end
                # no second stage but no early stopping or second stage but non-finite c or too low c
                if isfinite(c) & (n == n1)
                    @constraint(m, y[x1, n, c] == 0)
                end
                if !isfinite(c) & (n > n1)
                    @constraint(m, y[x1, n, c] == 0)
                end
                if !possible(n1, n, ss) # sample space constraints
                    @constraint(m,
                        y[x1, n, c] == 0
                    )
                end
            end
        end
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
        sum(dbinom(x1, n1, params.pess)*n*y[x1, n, c] for
            x1 in 0:n1,
            n  in nvals,
            c  in cvals
        )
    )
    return m, y
end

function _isfeasible(design::BinaryTwoStageDesign, params::MinimalMinimalExpectedSampleSize)
    all(power.(design, linspace(0, null(params))) .<= alpha(params) + .001) ? nothing : throw(InexactError())
    all(power.(design, linspace(alternative(params), 1)) .>= 1 - beta(params) - .001) ? nothing : throw(InexactError())
    return true
end
