function expected_sample_size(design::BinaryTwoStageDesign, prior::Function)::Float64
    @assert abs(quadgk(p -> prior(p), 0, 1, abstol = .001)[1] - 1.0) < .01
    return quadgk(p -> mean(SampleSize(design, p)*prior(p), 0, 1, abstol = .001)[1]
end

type MinimalExpectedSampleSizeVagueAlternative <: VagueAlternative
    n1range
    nmax
    alpha
    beta
    p0
    prior
    score
    solver
    function MinimalExpectedSampleSizeVagueAlternative{T1<:Integer, T2<:Real}(
        n1range::Union{Range{T1}, Vector{T1}},
        nmax::T1,
        alpha::T2,
        beta::T2,
        p0::T2,
        prior::Function,
        solver::MathProgBase.AbstractMathProgSolver
    )
        @assert all(n1range .>= 0)
        @assert all(n1range .<= nmax)
        @assert all([alpha; beta; p0] .>= 0.0)
        @assert all([alpha; beta; p0] .<= 1.0)
        new(
            n1range, nmax,
            alpha, beta,
            p0, prior,
            design -> expected_sample_size(design, prior)),
            solver
        )
    end
end

function _createProblem(
    n1::Int64,      # stage one sample size
    parameters::MinimalExpectedSampleSizeVagueAlternative
)
    # for brevity, extract planning parameters:
    nmax   = parameters.nmax
    alpha  = parameters.alpha
    beta   = parameters.beta
    p0     = parameters.p0
    prior  = parameters.prior
    # define base problem
    m, y, n1, nmax  = _createBaseProblem(n1, nmax) # c.f. util.jl
    # add type one error rate constraint
    @constraint(m,
        sum{
            dbinom(x1, n1, p0)*_cpr(x1, n1, n, c, p0)*y[x1, n, c],
            x1 = 0:n1,
            n = n1:nmax,
            c = [-Inf; 0:(nmax-1); Inf]
        } <= alpha
    )
    # conditional prior
    z::Float64 = quadgk(prior, p0, 1, .0001)[1]
    function cprior{T<:Real}(p::T)::Float64
        return prior(p)/z
    end
    # add type two error rate constraint (expected power)
    @constraint(m,
        sum{
            (quadgk(p -> dbinom(x1, n1, p)*cprior(p), p0, 1, abstol = 0.001)[1]*
                quadgk(p -> _cpr(x1, n1, n, c, p)*cprior(p), p0, 1, abstol = 0.001)[1]*
                y[x1, n, c]),
            x1 = 0:n1,
            n = n1:nmax,
            c = [-Inf; 0:(nmax - 1); Inf]
        } >= 1 - beta
    )
    # add optimality criterion
    @objective(m, Min,
        sum{
            quadgk(p -> dbinom(x1, n1, p)*prior(p), 0, 1, abstol = 0.001)[1]*n*y[x1, n, c],
            x1 = 0:n1,
            n  = n1:nmax,
            c  = [-Inf; 0:(nmax - 1); Inf]
        }
    )
    return m, y, n1, parameters
end
