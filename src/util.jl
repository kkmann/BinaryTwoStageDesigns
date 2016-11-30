"""
Creates a basline JuMP model
with only minimal constraints (contiguous stopping, unimodal n(x1)) and
functional constraint
"""
function _createBaseProblem(n1, n_max)
    m = Model()
    @variable(
        m,
        y[x1 = 0:n1, n = n1:n_max, c = [-Inf; 0:(n_max - 1); Inf]],
        Bin
    )
    # functional constraint: exactly one non-zero entry in (y) per x1
    for x1 in 0:n1
        @constraint(
            m,
            sum{y[x1, n, c], n = n1:n_max, c = [-Inf; 0:(n_max - 1); Inf]} == 1
        )
    end
    # several consistency requirements
    for x1 in 0:n1
        for n in n1:n_max
            for c in [-Inf; 0:(n_max - 1); Inf]
                # no second stage but no early stopping or second stage but non-finite c or too low c
                if (isfinite(c) & (n == n1)) | (!isfinite(c) & (n > n1)) | ((c < x1) & (n > n1))
                    @constraint(
                        m,
                        y[x1, n, c] == 0
                    )
                end
            end
        end
    end
    # implement contiguous regions for stopping for futility
    for x1 in 1:n1
        # if y[x1, n1, Inf] = 1 (stopping for futility at x1),
        # then y[x1 - 1, n1, Inf] = 1 must hold too (for x1 = 1:n1)
        # if y[x1, n1, Inf] = 0 the constraint is trivially fulfilled
        @constraint(
            m,
            y[x1 - 1, n1, Inf] - y[x1, n1, Inf] >= 0
        )
    end
    # implement contiguous regions for stopping for efficacy
    for x1 in 0:(n1 - 1)
        # if y[x1, n1, -Inf] = 1 (stopping for efficacy at x1),
        # then y[x1 + 1, n1, -Inf] = 1 must hold too (for x1 = 0:(n1 - 1))
        # if y[x1, n1, -Inf] = 0 the constraint is trivially fulfilled
        @constraint(
            m,
            y[x1 + 1, n1, -Inf] - y[x1, n1, -Inf] >= 0
        )
    end
    # implement unimodality of n(x1)
    @variable(
        m,
        _is_mode[x1 = 0:n1],
        Bin
    )
    # first: positive increment constraints (all increments before x1 are positive, if c is finite)
    for x1 in 0:n1
        for x1_ in 1:x1
            @constraint(
                m,
                sum{n*(y[x1_, n, c] - y[x1_ - 1, n, c]), n = (n1 + 1):(n_max), c = 0:(n_max - 1)} - 3*n_max*_is_mode[x1] >= -3*n_max
            )
        end
    end
    # second: negative increments (all increments after x1 are negative, if c is finite)
    for x1 in 0:n1
        for x1_ in (x1 + 1):(n1)
            @constraint(
                m,
                sum{n*(y[x1_, n, c] - y[x1_ - 1, n, c]), n = (n1 + 1):(n_max), c = 0:(n_max - 1)} + 3*n_max*_is_mode[x1] <= 3*n_max
            )
        end
    end
    @constraint(
        m,
        sum{_is_mode[x1], x1 = 0:n1} >= 1 # at least one mode (can be on boundary as well!)
    )
    return m, y, n1, n_max
end

function _extractSolution(y, n1, nmax)
    # extract solution
    n_func = zeros(Int64, n1 + 1)
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
    try
        design = BinaryTwoStageDesign(n_func, c_func)
    catch e
        println(e)
        error("could not create design from optimization result, consider changing solver or solver parameters to increase numerical accuracy.")
    end
end

function _cpr(x1, n1, n, c, p)
    # conditional probability to reject
    if x1 > c
        return(1.0)
    end
    if n - n1 + x1 <= c
        return(0.0)
    end
    return 1 - Distributions.cdf(Distributions.Binomial(n - n1, p), convert(Int64, c - x1))
end

function _addConditionalTypeOneErrorRateStageOneConstraint{T<:Parameters}(m, y, n1obs, x1obs, nna1, design, parameters::T)
    n1old = getInterimSampleSize(design)
    nmax  = parameters.nmax
    if n1old < n1obs
        for x1 in x1obs:(x1obs + nna1)
            @constraint(m,
                sum{
                    Distributions.pdf(Distributions.Binomial(n1obs - n1old, parameters.p0), j)*_cpr(x1 + j, n1obs, n, c, parameters.p0)*y[x1 + j, n, c],
                    j = 0:(n1obs - n1old),
                    n = n1obs:nmax,
                    c = [-Inf; 0:(nmax - 1); Inf]
                } <= _cpr(x1, n1old, getSampleSize(design, x1), getRejectionBoundary(design, x1), parameters.p0)
            )
        end
    end
    if n1old > n1obs
        for x1 in x1obs:(x1obs + nna1)
            @constraint(m,
                sum{
                    _cpr(x1, n1obs, n, c, parameters.p0)*y[x1, n, c],
                    n = n1obs:nmax,
                    c = [-Inf; 0:(nmax - 1); Inf]
                } <= sum([Distributions.pdf(Distributions.Binomial(n1old - n1obs, parameters.p0), j)*_cpr(x1 + j, n1old, getSampleSize(design, x1 + j), getRejectionBoundary(design, x1 + j), parameters.p0) for j in 0:(n1old - n1obs)])
            )
        end
    end
    return m, y
end
function _addConditionalTypeOneErrorRateStageTwoConstraint{T<:Parameters}(m, y, n1, x1obs, nna1, xobs, nna, design, parameters::T)
    n1old = getInterimSampleSize(design)
    nmax  = parameters.nmax
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
            return 1 - Distributions.cdf(Distributions.Binomial(n_x1 - design.n[x1 + 1], design.parameters.p0), c_x1 - x1 - l) # TODO: replace
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
            return 1 - Distributions.cdf(Distributions.Binomial(design.n[x1 + 1] - n_x1, design.parameters.p0), design.c[x1 + 1] - x1 - l) # TODO: replace
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
    return m, y
end
function _addConditionalTypeTwoErrorRateStageOneConstraint{T<:PointAlternative}(m, y, n1obs, x1obs, nna1, design, parameters::T)
    n1old = getInterimSampleSize(design)
    nmax  = parameters.nmax
    if n1old < n1obs
        for x1 in x1obs:(x1obs + nna1)
            @constraint(m,
                sum{
                    Distributions.pdf(Distributions.Binomial(n1obs - n1old, parameters.p1), j)*_cpr(x1 + j, n1obs, n, c, parameters.p1)*y[x1 + j, n, c],
                    j = 0:(n1obs - n1old),
                    n = n1obs:nmax,
                    c = [-Inf; 0:(nmax - 1); Inf]
                } >= _cpr(x1, n1old, getSampleSize(design, x1), getRejectionBoundary(design, x1), parameters.p1)
            )
        end
    end
    if n1old > n1obs
        for x1 in x1obs:(x1obs + nna1)
            @constraint(m,
                sum{
                    _cpr(x1, n1obs, n, c, parameters.p1)*y[x1, n, c],
                    n = n1obs:nmax,
                    c = [-Inf; 0:(nmax - 1); Inf]
                } >= sum([Distributions.pdf(Distributions.Binomial(n1old - n1obs, parameters.p1), j)*_cpr(x1 + j, n1old, getSampleSize(design, x1 + j), getRejectionBoundary(design, x1 + j), parameters.p1) for j in 0:(n1old - n1obs)])
            )
        end
    end
    return m, y
end

function _addInvarianceUnderImputationStageOneConstraint{T<:Parameters}(m, y, n1obs, x1obs, nna1, design, parameters::T)
    nmax  = parameters.nmax
    if nna1 > 0
        for x1 in (x1obs + 1):(x1obs + nna1)
            for n in n1obs:nmax
                @constraint(m,
                    sum{y[x1, n, c], c = [-Inf; 0:(nmax - 1); Inf]} - sum{y[x1 - 1, n, c], c = [-Inf; 0:(nmax - 1); Inf]} == 0
                )
            end
        end
    end
    return m, y
end

function _support(design::AbstractBinaryTwoStageDesign)
    n1     = getInterimSampleSize(design)
    nmax   = maximum(getSampleSize(design))
    return [[x1, x2] for x1 in 0:n1, x2 in 0:(nmax - n1) if x2 <= getSampleSize(opt, x1)] |> x-> hcat(x...)'
end

function _isInSupport{T<:Integer}(supp::Array{Int64, 2}, x1::T, x2::T)
    return any(map(i -> [x1; x2] == supp[i, :], 1:size(supp, 1)))
end
