type PointAlternative <: PlanningAssumptions
    n1range
    nmax
    alpha
    beta
    p0
    p1
    pess
    score
    solver
    function PointAlternative{T1<:Integer, T2<:Real}(
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
        new(n1range, nmax, alpha, beta, p0, p1, pess, design -> Distributions.mean(FinalSampleSize(design, pess)), solver)
    end
end

function _createProblem(
    n1::Int64,      # stage one sample size
    params::PointAlternative
)
    # for brevity, extract planning parameters:
    nmax   = params.nmax
    alpha  = params.alpha
    beta   = params.beta
    p0     = params.p0
    p1     = params.p1
    solver = params.solver
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
    return m, y, n1, params
end


function adaptToInterim{T<:Integer}(design::BinaryTwoStageDesign, n1obs::T, x1obs::T, nna1::T, params::PointAlternative) # isnt x1obs!
    # for brevity, extract planning parameters:
    nmax   = params.nmax
    alpha  = params.alpha
    beta   = params.beta
    p0     = params.p0
    p1     = params.p1
    solver = params.solver
    @assert n1obs <= maximum(params.n1range)
    # recreate the same problem conditional on the observed stage one sample size
    m, y, n1obs, params = _createProblem(n1obs, params)

    n1old = interimSampleSize(design)
    # case distinction for n1 <=> n1observed
    if n1old == n1obs
        # nothing to do
        return(design)
    end
    if n1old < n1obs
        for x1 in x1obs:(x1obs + nna1)
            @constraint(m,
                sum{
                    pdf(Binomial(n1obs - n1old, p0), j)*_cpr(x1 + j, n1obs, n, c, p0)*y[x1 + j, n, c],
                    j = 0:(n1obs - n1old),
                    n = n1obs:nmax,
                    c = [-Inf; 0:(nmax - 1); Inf]
                } <= _cpr(x1, n1old, design.n[x1 + 1], design.c[x1 + 1], p0)
            )
            @constraint(m,
                sum{
                    pdf(Binomial(n1obs - n1old, p1), j)*_cpr(x1 + j, n1obs, n, c, p1)*y[x1 + j, n, c],
                    j = 0:(n1obs - n1old),
                    n = n1obs:nmax,
                    c = [-Inf; 0:(nmax - 1); Inf]
                } >= _cpr(x1, n1old, design.n[x1 + 1], design.c[x1 + 1], p1)
            )
        end
    end
    if n1old > n1obs
        for x1 in x1obs:(x1obs + nna1)
            @constraint(m,
                sum{
                    _cpr(x1, n1obs, n, c, p0)*y[x1, n, c],
                    n = n1obs:nmax,
                    c = [-Inf; 0:(nmax - 1); Inf]
                } <= sum([pdf(Binomial(n1old - n1obs, p0), j)*_cpr(x1 + j, n1old, design.n[x1 + j + 1], design.c[x1 + j + 1], p0) for j in 0:(n1old - n1obs)])
            )
            println(x1)
            @constraint(m,
                sum{
                    _cpr(x1, n1obs, n, c, p1)*y[x1, n, c],
                    n = n1obs:nmax,
                    c = [-Inf; 0:(nmax - 1); Inf]
                } >= sum([pdf(Binomial(n1old - n1obs, p1), j)*_cpr(x1 + j, n1old, design.n[x1 + j + 1], design.c[x1 + j + 1], p1) for j in 0:(n1old - n1obs)])
            )
        end
    end
    # na constraints: x1obs .. x1obs + nna1 must have equal n!
    if nna1 > 0
        for x1 in (x1obs + 1):(x1obs + nna1)
            for n in n1obs:nmax
                @constraint(m,
                    sum{y[x1, n, c], c = [-Inf; 0:(nmax - 1); Inf]} - sum{y[x1 - 1, n, c], c = [-Inf; 0:(nmax - 1); Inf]} == 0
                )
            end
        end
    end
    setsolver(m, params.solver)
    status = solve(m)
    if !(status in (:Optimal, :UserLimit)) # no valid solution found!
        error("no feasible solution reached")
    end
    # extract solution
    design = _extractSolution(y, n1obs, params.nmax) # c.f. util.jl
    # TODO: check whether constraints also hold for p < p0!
    return design
end

function adaptToFinal{T<:Integer}(design::BinaryTwoStageDesign, x1obs::T, nna1::T, nobs::T, xobs::T, nna::T, params::PointAlternative)
    # for brevity, extract planning parameters:
    nmax   = params.nmax
    alpha  = params.alpha
    beta   = params.beta
    p0     = params.p0
    p1     = params.p1
    solver = params.solver
    @assert nobs <= nmax
    n1 = interimSampleSize(design)
    # recreate initial problem
    m, y, n1, params = _createProblem(n1, params)
    # case distinction for n <=> nobs
    if nobs == finalSampleSize(design)
        # nothing to do
        return(design)
    end
    # na constraints: x1obs .. x1obs + nna1 must have equal n!
    if nna1 > 0
        for x1 in (x1obs + 1):(x1obs + nna1)
            for n in n1:nmax
                @constraint(m,
                    sum{y[x1, n, c], c = [-Inf; 0:(nmax - 1); Inf]} - sum{y[x1 - 1, n, c], c = [-Inf; 0:(nmax - 1); Inf]} == 0
                )
            end
        end
    end
    # lhs
    function lhs(x1::Int64, l::Int64, n_x1::Int64, c_x1::Float64) # c must accept float for +- Inf
        if l > min(n_x1, design.n[x1 + 1]) - n1 # not possible
            return 0.0
        end
        if n_x1 <= design.n[x1 + 1] # new design has smaller value
            if x1 + l <= c_x1
                return 0.0
            else
                return 1.0
            end
        else
            # catch c = +- Inf
            if c_x1 == -Inf
                return(1.0)
            end
            if c_x1 == Inf
                return(0.0)
            end
            return 1 - Distributions.cdf(Distributions.Binomial(n_x1 - design.n[x1 + 1], p0), c_x1 - x1 - l) # TODO: replace
        end
    end
    # rhs
    function rhs(x1::Int64, l::Int64, n_x1::Int64, c_x1::Float64)
        if l > min(n_x1, design.n[x1 + 1]) - n1 # not possible
            return 0.0
        end
        if n_x1 < design.n[x1 + 1] # new design has smaller value
            # catch c = +- Inf
            if design.c[x1 + 1] == -Inf
                return(1.0)
            end
            if design.c[x1 + 1] == Inf
                return(0.0)
            end
            return 1 - Distributions.cdf(Distributions.Binomial(design.n[x1 + 1] - n_x1, p0), design.c[x1 + 1] - x1 - l) # TODO: replace
        else
            if x1 + l <= design.c[x1 + 1]
                return 0.0
            else
                return 1.0
            end
        end
    end
    for x1 in x1obs:(x1obs + nna1)
        for l in (xobs - x1obs):(xobs - x1obs + nna - nna1)
            @constraint(m,
                sum{
                    y[x1, n, c]*(lhs(x1, l, n, c) - rhs(x1, l, n, c)),
                    n = n1:nmax,
                    c = [-Inf; 0:(nmax - 1); Inf]
                } <= 0.0
            )
        end
    end
    # add constraints that ensures n(x1obs) = nobs (and for all other possible outcomes of stage one)
    for x1 in x1obs:(x1obs + nna1)
        @constraint(m,
            sum{y[x1, nobs, c], c = [-Inf; 0:(nmax -1); Inf]} >= 1
        )
    end
    setsolver(m, solver)
    status = solve(m)
    if !(status in (:Optimal, :UserLimit)) # no valid solution found!
        error("no feasible solution reached")
    end
    # get solution
    design = _extractSolution(y, n1, params.nmax)
    return design
end
