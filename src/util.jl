"""
Creates a basline JuMP model
with only minimal constraints (contiguous stopping, unimodal n(x1)) and
functional constraint
"""
function _createBaseProblem(n1, params)
    reg = regularization(params)
    eff = efficacy(params)
    ss = samplespace(params)
    !possible(n1, ss) ? throw(InexactError()) : nothing
    nmax = maxsamplesize(ss)
    nvals = n1:nmax
    if eff == StoppingForEfficacy
        cvals = [-Inf; 0:(nmax - 1); Inf]
    end
    if eff == NoStoppingForEfficacy
        cvals = [0:(nmax - 1); Inf]
    end
    m = Model()
    @variable(m, # indicator variables y[x1, n, c] == 1 iff n(x1) = n, c(x1) = c
        y[x1 = 0:n1, n = nvals, c = cvals],
        Bin
    )
    if reg == Unimodal
        @variable(m, # dummy variables for unimodality
            _is_mode[x1 = 0:n1],
            Bin
        )
    end
    if reg == GroupSequential
        @variable(m, # dummy variables for sample size
            samplesize[n = nvals],
            Bin
        )
        @constraint(m,
            sum(samplesize[n] for n in nvals) == 1
        )
        @variable(m, # dummy variables for critical value
            criticalvalue[c = cvals],
            Bin
        )
        @constraint(m,
            sum(criticalvalue[c] for c in 0:(nmax - 1)) == 1
        )
    end
    for x1 in 0:n1
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
        if (x1 < n1) & (eff == StoppingForEfficacy)
            # if y[x1, n1, -Inf] = 1 (stopping for eff at x1),
            # then y[x1 + 1, n1, -Inf] = 1 must hold too (for x1 = 0:(n1 - 1))
            # if y[x1, n1, -Inf] = 0 the constraint is trivially fulfilled
            @constraint(m,
                y[x1 + 1, n1, -Inf] - y[x1, n1, -Inf] >= 0
            )
        end
        for n in nvals
            for c in cvals
                # no second stage but no early stopping or second stage but non-finite c or too low c
                if (isfinite(c) & (n == n1)) | (!isfinite(c) & (n > n1)) | ((c < x1) & (n > n1))
                    @constraint(m,
                        y[x1, n, c] == 0
                    )
                end
                if !possible(n1, n, ss) # sample space constraints
                    @constraint(m,
                        y[x1, n, c] == 0
                    )
                end
                if (reg == GroupSequential) & isfinite(c)
                    @constraint(m, # n must equal the common sample size upon continuation
                        y[x1, n, c] - samplesize[n] == 0
                    )
                    @constraint(m, # c must equal the common critical value upon continuation
                        y[x1, n, c] - criticalvalue[c] == 0
                    )
                end
            end
        end
        if reg == Unimodal
            for x1_ in 1:x1
                @constraint(m,
                    sum(n*(y[x1_, n, c] - y[x1_ - 1, n, c]) for
                        n in (n1 + 1):nmax,
                        c in 0:(nmax - 1)
                    ) - 3*nmax*_is_mode[x1] >= -3*nmax
                )
            end
            for x1_ in (x1 + 1):n1
                @constraint(m,
                    sum(n*(y[x1_, n, c] - y[x1_ - 1, n, c]) for
                        n in (n1 + 1):nmax,
                        c in 0:(nmax - 1)
                    ) + 3*nmax*_is_mode[x1] <= 3*nmax
                )
            end
        end
    end
    if reg == Unimodal
        @constraint(m, # unimodality
            sum(_is_mode[x1] for x1 in 0:n1) >= 1 # at least one mode (can be on boundary as well!)
        )
    end
    return m, y
end

function _extractSolution(y, n1, params)
    ss = samplespace(params)
    eff = efficacy(params)
    nmax = maxsamplesize(ss)
    nvals = n1:nmax
    if eff == StoppingForEfficacy
        cvals = [-Inf; 0:(nmax - 1); Inf]
    end
    if eff == NoStoppingForEfficacy
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
        println(e)
        error("could not create design from optimization result, consider changing solver or solver parameters to increase numerical accuracy.")
    end
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
    return [[x1, x2] for x1 in 0:n1, x2 in 0:(nmax - n1) if (x2 <= getSampleSize(design, x1) - n1) & ((getSampleSize(design, x1) > n1) | (x2 == 0))] |> x-> hcat(x...)'
end

function _isInSupport{T<:Integer}(supp::Array{Int64, 2}, x1::T, x2::T)
    return any(map(i -> [x1; x2] == supp[i, :], 1:size(supp, 1)))
end

function dbinom{T1<:Integer, T2<:Real}(k::T1, n::T1, p::T2)::T2
    return Distributions.pdf(Distributions.Binomial(n, p), k)
end
