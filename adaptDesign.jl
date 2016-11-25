
function adaptToInterim(design:BinaryTwoStageDesign, n1_obs, x1_obs, nna1, params) # isnt x1_obs!
    # for brevity, extract planning parameters:
    nmax  = params.nmax
    alpha  = params.alpha
    beta   = params.beta
    p0    = params.p0
    p1    = params.p1
    solver = params.solver

    m, y, n1_obs, params = _create_point_alternative_problem(n1_obs, params)

    n1_old = length(design.n) - 1
    # case distinction for n1 <=> n1_observed
    if n1_old == n1_obs
        # nothing to do
        return(design)
    end
    if n1_old < n1_obs
        for x1 in x1_obs:(x1_obs + nna1)
            @constraint(m,
                sum{
                    pdf(Binomial(n1_obs - n1_old, p0), j)*_cpr(x1 + j, n1_obs, n, c, p0)*y[x1 + j, n, c],
                    j = 0:(n1_obs - n1_old),
                    n = n1_obs:nmax,
                    c = [-Inf; 0:(nmax - 1); Inf]
                } <= _cpr(x1, n1_old, design.n[x1 + 1], design.c[x1 + 1], p0)
            )
            @constraint(m,
                sum{
                    pdf(Binomial(n1_obs - n1_old, p1), j)*_cpr(x1 + j, n1_obs, n, c, p1)*y[x1 + j, n, c],
                    j = 0:(n1_obs - n1_old),
                    n = n1_obs:nmax,
                    c = [-Inf; 0:(nmax - 1); Inf]
                } >= _cpr(x1, n1_old, design.n[x1 + 1], design.c[x1 + 1], p1)
            )
        end
    end
    if n1_old > n1_obs
        for x1 in x1_obs:(x1_obs + nna1)
            @constraint(m,
                sum{
                    _cpr(x1, n1_obs, n, c, p0)*y[x1, n, c],
                    n = n1_obs:nmax,
                    c = [-Inf; 0:(nmax - 1); Inf]
                } <= sum([pdf(Binomial(n1_old - n1_obs, p0), j)*_cpr(x1 + j, n1_old, design.n[x1 + j + 1], design.c[x1 + j + 1], p0) for j in 0:(n1_old - n1_obs)])
            )
            println(x1)
            @constraint(m,
                sum{
                    _cpr(x1, n1_obs, n, c, p1)*y[x1, n, c],
                    n = n1_obs:nmax,
                    c = [-Inf; 0:(nmax - 1); Inf]
                } >= sum([pdf(Binomial(n1_old - n1_obs, p1), j)*_cpr(x1 + j, n1_old, design.n[x1 + j + 1], design.c[x1 + j + 1], p1) for j in 0:(n1_old - n1_obs)])
            )
        end
    end
    # na constraints: x1_obs .. x1_obs + nna1 must have equal n!
    if nna1 > 0
        for x1 in (x1_obs + 1):(x1_obs + nna1)
            for n in n1_obs:nmax
                @constraint(m,
                    sum{y[x1, n, c], c = [-Inf; 0:(nmax - 1); Inf]} - sum{y[x1 - 1, n, c], c = [-Inf; 0:(nmax - 1); Inf]} == 0
                )
            end
        end
    end
    setsolver(m, solver)
    status = solve(m)
    if !(status in (:Optimal, :UserLimit)) # no valid solution found!
        error("no feasible solution reached")
    end
    # extract solution
    n_func = zeros(Float64, n1_obs + 1)
    c_func = zeros(Float64, n1_obs + 1) # need float for +/- Inf
    for x1 in 0:n1_obs
        for n in n1_obs:nmax
            for c in [-Inf; 0:(nmax - 1); Inf]
                if getvalue(y[x1, n, c]) == 1
                    n_func[x1 + 1] = n
                    c_func[x1 + 1] = c
                end
            end
        end
    end
    # check whether constraints also hold for p < p0!
    design = BinaryTwoStageDesign(n_func, c_func)
    return design
end

function adaptToFinal(design::BinaryTwoStageDesign, x1_obs, nna1, n_obs, x_obs, nna, params)
    # for brevity, extract planning parameters:
    nmax  = params.nmax
    alpha  = params.alpha
    beta   = params.beta
    p0    = params.p0
    p1    = params.p1
    solver = params.solver

    if n_obs > nmax
        error("n_obs must be smaller or equal to nmax!")
    end

    n1 = length(design.n) - 1

    m, y, n1, params = _create_point_alternative_problem(n1, params)

    # case distinction for n <=> n_obs
    if n_obs == design.n[x1_obs + 1]
        # nothing to do
        return(design)
    end
    # na constraints: x1_obs .. x1_obs + nna1 must have equal n!
    if nna1 > 0
        for x1 in (x1_obs + 1):(x1_obs + nna1)
            for n in n1:nmax
                @constraint(m,
                    sum{y[x1, n, c], c = [-Inf; 0:(nmax - 1); Inf]} - sum{y[x1 - 1, n, c], c = [-Inf; 0:(nmax - 1); Inf]} == 0
                )
            end
        end
    end
    println("hi")
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
            return 1 - cdf(Binomial(n_x1 - design.n[x1 + 1], p0), c_x1 - x1 - l)
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
            return 1 - cdf(Binomial(design.n[x1 + 1] - n_x1, p0), design.c[x1 + 1] - x1 - l)
        else
            if x1 + l <= design.c[x1 + 1]
                return 0.0
            else
                return 1.0
            end
        end
    end
    println("hey")
    for x1 in x1_obs:(x1_obs + nna1)
        for l in (x_obs - x1_obs):(x_obs - x1_obs + nna - nna1)
            @constraint(m,
                sum{
                    y[x1, n, c]*(lhs(x1, l, n, c) - rhs(x1, l, n, c)),
                    n = n1:nmax,
                    c = [-Inf; 0:(nmax - 1); Inf]
                } <= 0.0
            )
        end
    end
    println("ho")
    # add constraint that ensures n(x1_obs) = n_obs (and for all other possible outcomes of stage one)
    @constraint(m,
        sum{y[x1_obs, n_obs, c], c = [-Inf; 0:(nmax -1); Inf]} >= 1
    )
    setsolver(m, solver)
    status = solve(m)
    if !(status in (:Optimal, :UserLimit)) # no valid solution found!
        error("no feasible solution reached")
    end
    # extract solution
    n_func = zeros(Float64, n1 + 1)
    c_func = zeros(Float64, n1 + 1) # need float for +/- Inf
    for x1 in 0:n1
        for n in n1:nmax
            for c in [-Inf; 0:(nmax - 1); Inf]
                if getvalue(y[x1, n, c]) == 1
                    n_func[x1 + 1] = n
                    c_func[x1 + 1] = c
                end
            end
        end
    end
    # check whether constraints also hold for p < p0!
    design = BinaryTwoStageDesign(n_func, c_func)
    return design
end
