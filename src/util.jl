"""
Any planning assumption must at least contain the fields n1range as the
allowable range of n1 values, nmax as maximal final n and score(design) as
function giving the score of a design (smaller = better).
"""
abstract PlanningAssumptions

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
    return BinaryTwoStageDesign(n_func, c_func)
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
