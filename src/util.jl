"""
Creates a basline JuMP model
with only minimal constraints (contiguous stopping, unimodal n(x1)) and
functional constraint TODO make samplespace restriction more effective by using maxsamplesize(ss, n1)
"""
function _createBaseProblem(n1, params) # regularize by penalizing total variation  of c and n
    ss  = samplespace(params)
    !possible(n1, ss) ? throw(InexactError()) : nothing
    nmax = maxsamplesize(ss)
    nvals = n1:nmax
    cvalsfinite = 0:(nmax - 1)
    if allowsstoppingforefficacy(params)
        cvals = [-Inf; cvalsfinite; Inf]
        cvalsinfinite = [-Inf; Inf]
    else
        cvals = [cvalsfinite; Inf]
        cvalsinfinite = [Inf]
    end

    m = Model()
    @variable(m, # indicator variables y[x1, n, c] == 1 iff n(x1) = n, c(x1) = c
        y[x1 = 0:n1, n = nvals, c = cvals],
        Bin
    )
    if !isgroupsequential(params)
        @variable(m, # dummy variables for unimodality
            _is_mode[x1 = 0:n1],
            Bin
        )
    else
        @variable(m, # dummy variables for sample size and critical value
            cont[x1 = 0:n1],
            Bin
        )
    end
    for x1 in 0:n1
        if isgroupsequential(params)
            @constraint(m, # lhs is one if design continues at x1
                sum(y[x1, n, c] for n in nvals, c in cvalsinfinite) - (1 - cont[x1]) == 0
            )
        end
        @constraint(m, # functional constraint: exactly one non-zero entry in (y) per x1
            sum(y[x1, n, c] for n in nvals, c in cvals) == 1
        )
        # if y[x1, n1, Inf] = 1 (stopping for futility at x1),
        # then y[x1 - 1, n1, Inf] = 1 must hold too (for x1 = 1:n1)
        # if y[x1, n1, Inf] = 0 the constraint is trivially fulfilled
        if x1 > 0
            @constraint(m,
                y[x1 - 1, n1, Inf] - y[x1, n1, Inf] >= 0
            )
        end
        if (x1 < n1) & allowsstoppingforefficacy(params)
            # if y[x1, n1, -Inf] = 1 (stopping for eff at x1),
            # then y[x1 + 1, n1, -Inf] = 1 must hold too (for x1 = 0:(n1 - 1))
            # if y[x1, n1, -Inf] = 0 the constraint is trivially fulfilled
            @constraint(m,
                y[x1 + 1, n1, -Inf] - y[x1, n1, -Inf] >= 0
            )
        end
        for n in nvals
            for c in cvals
                if isfinite(c)
                    if c >= n
                        @constraint(m, y[x1, n, c] == 0)
                    end
                    if !(isgroupsequential(params) & !allowsstoppingforefficacy(params))
                        # we need this restriction as group sequential with
                        # no stopping for efficacy otheriwse implies c >= n1
                        # which is senseless
                        if c - x1 >= n - n1 # x1 + x2 > c impossible
                            @constraint(m, y[x1, n, c] == 0)
                        end
                    end
                end
                # no second stage but no early stopping or second stage but non-finite c or too low c
                if isfinite(c) & (n == n1)
                    @constraint(m, y[x1, n, c] == 0)
                end
                # if !isfinite(c) & (n > n1) # not necessary if design may continue with decision after stage one
                #     @constraint(m, y[x1, n, c] == 0)
                # end
                if ( (c < x1) & (n > n1) ) & !(isgroupsequential(params) & !allowsstoppingforefficacy(params))
                    @constraint(m, y[x1, n, c] == 0)
                end
                if !possible(n1, n, c, ss) # sample space constraints
                    @constraint(m,
                        y[x1, n, c] == 0
                    )
                end
                if isgroupsequential(params) & isfinite(c)
                    if x1 < n1
                        @constraint(m, # n, c must be equal for all x1 with finite c
                            y[x1, n, c] - y[x1 + 1, n, c] - (2 - cont[x1] - cont[x1 + 1]) <= 0
                        )
                        @constraint(m, # n, c must be equal for all x1 with finite c
                            y[x1, n, c] - y[x1 + 1, n, c] + (2 - cont[x1] - cont[x1 + 1]) >= 0
                        )
                    end
                end
            end
        end
        if !isgroupsequential(params)
            for x1_ in 1:x1
                @constraint(m,
                    sum(n*(y[x1_, n, c] - y[x1_ - 1, n, c]) for
                        n in (n1 + 1):nmax,
                        c in cvalsfinite
                    ) - 3*nmax*_is_mode[x1] >= -3*nmax
                )
            end
            for x1_ in (x1 + 1):n1
                @constraint(m,
                    sum(n*(y[x1_, n, c] - y[x1_ - 1, n, c]) for
                        n in (n1 + 1):nmax,
                        c in cvalsfinite
                    ) + 3*nmax*_is_mode[x1] <= 3*nmax
                )
            end
        end
    end
    if !isgroupsequential(params)
        @constraint(m, # unimodality
            sum(_is_mode[x1] for x1 in 0:n1) >= 1 # at least one mode (can be on boundary as well!)
        )
    end
    return m, y
end

function _extractSolution(y, n1, params)
    ss = samplespace(params)
    nmax = maxsamplesize(ss)
    nvals = n1:nmax
    if allowsstoppingforefficacy(params)
        cvals = [-Inf; 0:(nmax - 1); Inf]
    else
        cvals = [0:(nmax - 1); Inf]
    end
    nvec = zeros(Int64, n1 + 1)
    cvec = zeros(Float64, n1 + 1) # need float for +/- Inf
    for x1 in 0:n1
        for n in nvals
            for c in cvals
                if getvalue(y[x1, n, c]) == 1
                    nvec[x1 + 1] = n
                    cvec[x1 + 1] = c
                end
            end
        end
    end
    try
        design = BinaryTwoStageDesign(nvec, cvec, params)
    catch e
        println(nvec)
        println(cvec)
        error("could not create design from optimization result, consider changing solver or solver parameters to increase numerical accuracy.")
    end
end


function _addconditionaltypeonestageoneconstraint(
    m, y, n1obs, x1obs, nna1, design::BinaryTwoStageDesign
)
    params = parameters(design)
    p0     = null(params)
    n1old  = interimsamplesize(design)
    nmax   = maxsamplesize(params)
    cvalsfinite = 0:(nmax - 1)
    if allowsstoppingforefficacy(params)
        cvals = [-Inf; cvalsfinite; Inf]
        cvalsinfinite = [-Inf; Inf]
    else
        cvals = [cvalsfinite; Inf]
        cvalsinfinite = [Inf]
    end
    if n1old < n1obs
        for x1 in x1obs:(x1obs + nna1)
            @constraint(m,
                sum(
                    dbinom(j, n1obs - n1old, p0)*_cpr(x1 + j, n1obs, n, c, p0)*y[x1 + j, n, c] for
                    j = 0:(n1obs - n1old),
                    n = n1obs:nmax,
                    c = cvals
                ) <= _cpr(x1, n1old, samplesize(design, x1), criticalvalue(design, x1), p0)
            )
        end
    end
    if n1old > n1obs
        for x1 in x1obs:(x1obs + nna1)
            @constraint(m,
                sum(
                    _cpr(x1, n1obs, n, c, null(params))*y[x1, n, c] for
                    n = n1obs:nmax,
                    c = [-Inf; 0:(nmax - 1); Inf]
                ) <= sum([dbinom(j, n1old - n1obs, p0)*_cpr(x1 + j, n1old, samplesize(design, x1 + j), criticalvalue(design, x1 + j), p0) for
                    j in 0:(n1old - n1obs)])
            )
        end
    end
    return m, y
end
function _addconditionaltypeonestagetwoconstraint(
    m, y, n1, x1obs, nna1, xobs, nna, design
)
    params = parameters(design)
    p0     = null(params)
    n1old  = interimsamplesize(design)
    nmax   = maxsamplesize(params)
    nvals  = n1old:nmax
    cvalsfinite = 0:(nmax - 1)
    if allowsstoppingforefficacy(params)
        cvals = [-Inf; cvalsfinite; Inf]
        cvalsinfinite = [-Inf; Inf]
    else
        cvals = [cvalsfinite; Inf]
        cvalsinfinite = [Inf]
    end
    function lhs(x1::Int64, l::Int64, n_x1::Int64, c_x1::Float64) # c must accept float for +- Inf
        if l > min(n_x1, samplesize(design, x1)) - n1 # not possible
            return 0.0
        end
        if n_x1 <= samplesize(design, x1) # new design has smaller value
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
            return 1 - Distributions.cdf(Distributions.Binomial(n_x1 - samplesize(design, x1), p0), c_x1 - x1 - l) # TODO: replace
        end
    end
    # rhs
    function rhs(x1::Int64, l::Int64, n_x1::Int64, c_x1::Float64)
        if l > min(n_x1, samplesize(design, x1)) - n1 # not possible
            return 0.0
        end
        if n_x1 < samplesize(design, x1) # new design has smaller value
            # catch c = +- Inf
            if samplesize(design, x1) == -Inf
                return(1.0)
            end
            if criticalvalue(design, x1) == Inf
                return(0.0)
            end
            return 1 - Distributions.cdf(Distributions.Binomial(samplesize(design, x1) - n_x1, p0), samplesize(design, x1) - x1 - l) # TODO: replace
        else
            if x1 + l <= samplesize(design, x1)
                return 0.0
            else
                return 1.0
            end
        end
    end
    for x1 in x1obs:(x1obs + nna1)
        for l in (xobs - x1obs):(xobs - x1obs + nna - nna1)
            @constraint(m,
                sum(y[x1, n, c]*(lhs(x1, l, n, c) - rhs(x1, l, n, c)) for
                    n in nvals,
                    c in cvals
                ) <= 0.0
            )
        end
    end
    return m, y
end

function _addinvarianceimputationstageoneconstraint(
    m, y, n1obs, x1obs, nna1, design
)
    params = parameters(design)
    nmax   = maxsamplesize(params)
    cvalsfinite = 0:(nmax - 1)
    if allowsstoppingforefficacy(params)
        cvals = [-Inf; cvalsfinite; Inf]
        cvalsinfinite = [-Inf; Inf]
    else
        cvals = [cvalsfinite; Inf]
        cvalsinfinite = [Inf]
    end
    if nna1 > 0
        for x1 in (x1obs + 1):(x1obs + nna1)
            for n in n1obs:nmax
                @constraint(m,
                    sum(y[x1, n, c] for
                        c = cvals
                    ) - sum(y[x1 - 1, n, c] for
                        c = cvals
                    ) == 0
                )
            end
        end
    end
    return m, y
end


#
# function _isInSupport{T<:Integer}(supp::Array{Int64, 2}, x1::T, x2::T)
#     return any(map(i -> [x1; x2] == supp[i, :], 1:size(supp, 1)))
# end

dbinom(k, n, p) = Distributions.pdf(Distributions.Binomial(n, p), k)

qnorm(p) = Distributions.quantile(Distributions.Normal(0, 1), p)
