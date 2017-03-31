"""
Creates a basline JuMP model
with only minimal constraints (contiguous stopping, unimodal n(x1)) and
functional constraint
"""
function _createBaseProblem(n1, params) 
    ss  = samplespace(params)
    !possible(n1, ss) ? throw(InexactError()) : nothing
    nmax = maxsamplesize(ss, n1)
    nvals = getnvals(ss, n1)
    cvalsfinite = 0:(maximum(nvals) - 1)
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
        # contiguous stopping for futility
        if x1 > 0
            @constraint(m,
                y[x1 - 1, n1, Inf] - y[x1, n1, Inf] >= 0
            )
        end
        if (x1 < n1) & allowsstoppingforefficacy(params)
            # contiguous stopping for efficacy
            @constraint(m,
                y[x1 + 1, n1, -Inf] - y[x1, n1, -Inf] >= 0
            )
        end
        for n in nvals
            if isgroupsequential(params)
                @constraint(m, # groupsequential designs can only have fixed decision with early stopping
                    sum(y[x1, n1, c] for c in cvalsinfinite) == sum(y[x1, n, c] for n in nvals, c in cvalsinfinite)
                )
            end
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
                # c finite => n > n1
                if isfinite(c) & (n == n1)
                    @constraint(m, y[x1, n, c] == 0)
                end
                if (x1 > c) & (c != -Inf) & (!isgroupsequential(params) | (isgroupsequential(params) & allowsstoppingforefficacy(params))) # decision is fixed but c is not valid
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
                        n in nvals[2:end],
                        c in cvals
                    ) - 3*nmax*_is_mode[x1] >= -3*nmax
                )
            end
            for x1_ in (x1 + 1):n1
                @constraint(m,
                    sum(n*(y[x1_, n, c] - y[x1_ - 1, n, c]) for
                        n in nvals[2:end],
                        c in cvals
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
    ss    = samplespace(params)
    nmax  = maxsamplesize(ss, n1)
    nvals = getnvals(ss, n1)
    cvalsfinite = 0:(maximum(nvals) - 1)
    if allowsstoppingforefficacy(params)
        cvals = [-Inf; cvalsfinite; Inf]
        cvalsinfinite = [-Inf; Inf]
    else
        cvals = [cvalsfinite; Inf]
        cvalsinfinite = [Inf]
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
        return design
    catch e
        println("trying to recover solution") # sometimes values are set to 0 due to numerical issues
        println(nvec)
        println(cvec)
        try
            for i in 1:length(nvec)
                if nvec[i] == 0
                    leftn  = NaN
                    rightn = NaN
                    leftc  = NaN
                    rightc = NaN
                    if i > 1
                        j = 1
                        while (nvec[i - 1] == 0) & (i - j - 1 >= 1)
                            j -= 1
                            if nvec[i - j] != 0 # found left interpolation value
                                break
                            end
                        end
                        if nvec[i - j] != 0
                            leftn = nvec[i - j]
                            leftc = cvec[i - j]
                        end
                    end
                    if i < length(nvec)
                        j = 1
                        while (nvec[i + 1] == 0) & (i + j + 1 <= length(nvec))
                            j += 1
                            if nvec[i + j] != 0 # found right interpolation value
                                break
                            end
                        end
                        if nvec[i + j] != 0
                            rightn = nvec[i + j]
                            rightc = cvec[i + j]
                        end
                    end
                    if isnan(leftn) & isnan(rightn)
                        continue
                    end
                    if isnan(leftn)
                        nvec[i] = rightn
                        cvec[i] = rightc
                    end
                    if isnan(rightn)
                        nvec[i] = leftn
                        cvec[i] = leftc + 1 # need to account for x1 + 1
                    end
                    if !isnan(leftn) & !isnan(rightn)
                        nvec[i] = round((leftn + rightn)/2)
                        cvec[i] = max(leftc + 1, ceil((leftc + rightc)/2))
                    end
                end
            end
            println(nvec)
            println(cvec)
            design = BinaryTwoStageDesign(nvec, cvec, params)
            return design
        catch e
            println(e)
            error("could not create design from optimization result, consider changing solver or solver parameters to increase numerical accuracy.")
        end
    end
end


function _addconditionaltypeonestageoneconstraint(
    m, y, n1obs, x1obs, nna1, design::BinaryTwoStageDesign
)
    params = parameters(design)
    ss     = samplespace(params)
    nvals  = getnvals(ss, n1)
    if ss.stepsize != 1
        error("stepsize must be one")
    end
    p0     = null(params)
    n1old  = interimsamplesize(design)
    nmax   = maxsamplesize(params, n1)
    cvalsfinite = 0:(maximum(nvals) - 1)
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
    ss     = samplespace(params)
    nvals  = getnvals(ss, n1)
    if ss.stepsize != 1
        error("stepsize must be one")
    end
    p0     = null(params)
    n1old  = interimsamplesize(design)
    nmax   = maxsamplesize(params, n1)
    nvals  = n1old:nmax
    cvalsfinite = 0:(maximum(nvals) - 1)
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
    ss     = samplespace(params)
    nvals  = getnvals(ss, n1)
    if ss.stepsize != 1
        error("stepsize must be one")
    end
    nmax   = maxsamplesize(params, n1)
    cvalsfinite = 0:(maximum(nvals) - 1)
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

dbinom(k, n, p) = Distributions.pdf(Distributions.Binomial(n, p), k)

qnorm(p) = Distributions.quantile(Distributions.Normal(0, 1), p)

function findgrid(prior, l, u, n; resolution = 10000) # todo maybe put pivots in the middle of the intervals?
    candidates = collect(linspace(l, u, resolution + 2))[2:(resolution + 1)]
    quantiles  = collect(linspace(l, u, n + 2))[2:(n + 1)]
    modprior(p) = .9*prior(p) +  .1*Distributions.pdf(Distributions.Beta(.5, .5), p)
    cdfrelaxed = cumsum(modprior.(candidates))
    cdfrelaxed = cdfrelaxed / cdfrelaxed[resolution]
    cdf = cumsum(prior.(candidates))
    cdf = cdf / cdf[resolution]
    pivots = zeros(n)
    cdfpiv = zeros(n)
    i = 1
    for j in 1:n
        while (cdfrelaxed[i] < quantiles[j]) & (i < resolution)
            i += 1
        end
        pivots[j] = candidates[i]
        cdfpiv[j] = cdf[i]
    end
    pivots = (pivots .+ [l; pivots[1:(n - 1)]] ) / 2
    dcdf = cdfpiv - [0; cdfpiv[1:(n - 1)]] # vector of first order differences
    dcdf = dcdf/sum(dcdf)
    return pivots, dcdf
end
