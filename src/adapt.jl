function adapt(design::BinaryTwoStageDesign, data::DataFrames.DataFrame)
    # check data frame
    !isa(data[1], DataArrays.DataArray{Bool, 1}) ? error("first column of 'data' must be binary DataArray, true=response") : nothing
    !isa(data[2, 1], Integer) ? error("second column of 'data' must be integer DataArray, 1 = stage 1, 2 = stage 2") : nothing
    !all(map(x -> x in (1, 2), data[2])) ? error("second column of 'data' can only hold integers 1 or 2 (stages)") : nothing
    n1plan = interimsamplesize(design)
    n1obs::Integer  = data[1][data[2] .== 1] |> x -> x[!DataArrays.isna(x)] |> sum
    nna1::Integer   = data[1][data[2] .== 1] |> x -> x[DataArrays.isna(x)]  |> length
    n2obs::Integer  = data[1][data[2] .== 2] |> x -> x[!DataArrays.isna(x)] |> sum
    nna2::Integer   = data[1][data[2] .== 2] |> x -> x[DataArrays.isna(x)]  |> length
    nobs::Integer   = size(data, 1)
    if nobs == n1obs # only stage 1
        if (n1obs == n1plan) & (nna1 == 0) # nothing to do
            return design
        else
            !(n1obs in design.params.samplespace.n1range) ? error("observed n1 not compatible with sample space, consider modifying it") : nothing
            push!(design.params.samplespace.snvals, n1obs) # guarantee that observed n1 is in sample space
            push!(design.params.samplespace.snvals, nobs) # and observed sample size
            ipm = IPModel(samplespace(parameters(design)), n1obs)
            completemodel(imp, parameters(design), n1obs)
            addconditionaltypeoneerrorrateconstraint!(ipm, n1obs, x1obs, nna1, design)
            addinvarianceimputationconstraint!(ipm, n1obs, x1obs, nna1)
        end
    end
    if nobs > n1obs # also stage two (cannot be smaller!)
        if (n1obs == n1plan) & (nna1 == 0) & (nobs = samplesize(design, x1obs)) # we do not check for nna2 == 0 as this does not require adaptation, nothing to do!
            return design
        end
        addconditionaltypeoneerrorrateconstraint_stagetwo!(ipm, n1, x1obs, nna1, x1obs + x2obs, nna1 + nna2, design)
        for x1 in x1obs:(x1obs + nna1) # for all possible stage 1 outcome n must equal observed value
            @constraint(ipm.m,
                sum(ipm.y[x1, n1obs + n2obs, c] for c in cvals) >= 1
            )
        end
    end
    setsolver(ipm, solver)
    status = solve(ipm.m)
    if !(status in (:Optimal, :UserLimit)) # no valid solution found!
        error("no feasible solution reached")
    end
    newDesign = extractsolution(ipm, params)
    # # TODO: check whether constraints also hold for p < p0!
    return newDesign
end

function addconditionaltypeoneerrorrateconstraint!(
    ipm::IPModel, n1obs::Integer, x1obs::Integer, nna1::Integer, design::BinaryTwoStageDesign
)
    params = parameters(design)
    ss     = samplespace(params)
    n1old  = interimsamplesize(design)
    p0     = design |> parameters |> null
    nmax   = maxsamplesize(ss, n1obs)
    nvals  = getnvals(ss, n1obs)
    cvals =  [-Inf; getcvals(ss, n1obs); Inf]
    if n1old < n1obs
        for x1 in x1obs:(x1obs + nna1)
            @constraint(ipm.m,
                sum(
                    dbinom(j, n1obs - n1old, p0)*_cpr(x1 + j, n1obs, n, c, p0)*ipm.y[x1 + j, n, c] for
                    j = 0:(n1obs - n1old),
                    n = n1obs:nmax,
                    c = cvals
                ) <= _cpr(x1, n1old, samplesize(design, x1), criticalvalue(design, x1), p0)
            )
        end
    end
    if n1old > n1obs
        for x1 in x1obs:(x1obs + nna1)
            @constraint(ipm.m,
                sum(
                    _cpr(x1, n1obs, n, c, null(params))*ipm.y[x1, n, c] for
                    n = n1obs:nmax,
                    c = cvals
                ) <= sum([dbinom(j, n1old - n1obs, p0)*_cpr(x1 + j, n1old, samplesize(design, x1 + j), criticalvalue(design, x1 + j), p0) for
                    j in 0:(n1old - n1obs)])
            )
        end
    end
    return ipm
end
function addconditionaltypeoneerrorrateconstraint_stagetwo!(
    ipm, n1, x1obs, nna1, xobs, nna, design
)
    params = parameters(design)
    ss     = samplespace(params)
    n1old  = interimsamplesize(design)
    p0     = design |> parameters |> null
    nmax   = maxsamplesize(ss, n1obs)
    nvals  = getnvals(ss, n1obs)
    cvals =  [-Inf; getcvals(ss, n1obs); Inf]
    function lhs(x1, l, n_x1, c_x1) # c must accept float for +- Inf
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
    function rhs(x1, l, n_x1, c_x1)
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
            @constraint(ipm.m,
                sum(ipm.y[x1, n, c]*(lhs(x1, l, n, c) - rhs(x1, l, n, c)) for
                    n in nvals,
                    c in cvals
                ) <= 0.0
            )
        end
    end
    return ipm
end

function addinvarianceimputationconstraint!(
    ipm::IPModel, n1obs::Integer, x1obs::Integer, nna1::Integer
)
    nvals  = getnvals(ipm, n1obs)
    cvals =  [-Inf; getcvals(ipm, n1obs); Inf]
    if nna1 > 0
        for x1 in (x1obs + 1):(x1obs + nna1)
            for n in nvals
                @constraint(ipm.m,
                    sum(ipm.y[x1, n, c] for
                        c = cvals
                    ) - sum(y[x1 - 1, n, c] for
                        c = cvals
                    ) == 0
                )
            end
        end
    end
    return ipm
end

function adapt{T<:Integer, TS<:MathProgBase.AbstractMathProgSolver}(
    design::BinaryTwoStageDesign,
    outcome::DataArrays.DataVector{Bool},
    stage::Vector{T},
    solver::TS
)
    length(stage) == length(outcome) ? nothing : throw(InexactError())
    # check that parameters is not NoParameters (no optimization routine!)
    params = parameters(design)
    typeof(params) != NoParameters ? nothing : throw(InexactError())
    # check that stage is only 1, 2 (cannot be missing)
    all([s in [1 2] for s in stage]) ? nothing : throw(InexactError())
    length(stage) == length(outcome) ? nothing : throw(InexactError())
    # check that parameters is not NoParameters (no optimization routine!)
    params = parameters(design)
    typeof(params) != NoParameters ? nothing : throw(InexactError())
    # check that stage is only 1, 2 (cannot be missing)
    all([s in [1 2] for s in stage]) ? nothing : throw(InexactError())
    n1planned = interimsamplesize(design)
    nmax      = maxsamplesize(params)
    n1obs     = length(stage[stage .== 1])
    n1min     = min(n1planned, n1obs)
    x1obs     = outcome[stage .== 1] |> x -> x[!DataArrays.isna(x)] |> sum
    nna1      = sum(DataArrays.isna(outcome[stage .== 1][1:n1min]))
    xobs      = outcome |> x -> x[!DataArrays.isna(x)] |> sum
    nna       = sum(DataArrays.isna(outcome))
    nobs      = length(stage)
    if (n1planned == n1obs) & (nobs == samplesize(design, x1obs)) & (nna == 0) # nothing to do
        return design
    end
    # recreate the same problem conditional on the observed stage one sample size
    m, y = _createProblem(n1obs, params)
    m, y = _addconditionaltypeonestageoneconstraint(m, y, n1obs, x1obs, nna1, design)
    m, y = _addinvarianceimputationstageoneconstraint(m, y, n1obs, x1obs, nna1, design)
    if any(stage .== 2)
        m, y = _addconditionaltypeonestagetwoconstraint(m, y, n1planned, x1obs, nna1, xobs, nna, design)
        cvalsfinite = 0:(nmax - 1)
        if allowsstoppingforefficacy(params)
            cvals = [-Inf; cvalsfinite; Inf]
        else
            cvals = [cvalsfinite; Inf]
        end
        println(x1obs:(x1obs + nna1))
        for x1 in x1obs:(x1obs + nna1)
            @constraint(m,
                sum(y[x1, nobs, c] for c in cvals) >= 1
            )
        end
    end
    setsolver(m, solver)
    status = solve(m)
    if !(status in (:Optimal, :UserLimit)) # no valid solution found!
        error("no feasible solution reached")
    end
    # extract solution
    newDesign = _extractSolution(y, n1obs, params) # c.f. util.jl
    # # TODO: check whether constraints also hold for p < p0!
    return newDesign
end
