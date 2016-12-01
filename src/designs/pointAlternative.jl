type MinimalExpectedSampleSizePointAlternative <: PointAlternative
    n1range
    nmax
    alpha
    beta
    p0
    p1
    pess
    score
    solver
    function MinimalExpectedSampleSizePointAlternative{T1<:Integer, T2<:Real}(

        n1range::Union{Range{T1}, Vector{T1}},
        nmax::T1,
        alpha::T2,
        beta::T2,
        p0::T2,
        p1::T2,
        pess::T2,
        solver::MathProgBase.AbstractMathProgSolver
    )
        @assert all(n1range .>= 0)
        @assert all(n1range .<= nmax)
        @assert all([alpha; beta; p0; p1; pess] .>= 0.0)
        @assert all([alpha; beta; p0; p1; pess] .<= 1.0)
        @assert p0 < p1
        new(n1range, nmax, alpha, beta, p0, p1, pess, design -> mean(SampleSize(design, pess)), solver)
    end
end

function _createProblem(
    n1::Int64,      # stage one sample size
    parameters::MinimalExpectedSampleSizePointAlternative
)
    # for brevity, extract planning parameters:
    nmax   = parameters.nmax
    alpha  = parameters.alpha
    beta   = parameters.beta
    p0     = parameters.p0
    p1     = parameters.p1
    solver = parameters.solver
    # define base problem
    m, y, n1, nmax  = _createBaseProblem(n1, nmax) # c.f. util.jl
    # add type one error rate constraint
    @constraint(m,
        sum{
            Distributions.pdf(Distributions.Binomial(n1, p0), x1)*_cpr(x1, n1, n, c, p0)*y[x1, n, c],
            x1 = 0:n1,
            n = n1:nmax,
            c = [-Inf; 0:(nmax-1); Inf]
        } <= alpha
    )
    # add type two error rate constraint (power)
    @constraint(m,
        sum{
            Distributions.pdf(Distributions.Binomial(n1, p1), x1)*_cpr(x1, n1, n, c, p1)*y[x1, n, c],
            x1 = 0:n1,
            n = n1:nmax,
            c = [-Inf; 0:(nmax - 1); Inf]
        } >= 1 - beta
    )
    # add optimality criterion
    @objective(m, Min,
        sum{
            Distributions.pdf(Distributions.Binomial(n1, p1), x1)*n*y[x1, n, c],
            x1 = 0:n1,
            n  = n1:nmax,
            c  = [-Inf; 0:(nmax - 1); Inf]
        }
    )
    return m, y, n1, parameters
end
