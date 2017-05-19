function adapt(design::BinaryTwoStageDesign, data::DataFrames.DataFrame, solver; eta = 10000)
    # check data frame
    !isa(data[1], DataArrays.DataArray{Bool, 1}) ? error("first column of 'data' must be binary DataArray, true=response") : nothing
    !isa(data[2, 1], Integer) ? error("second column of 'data' must be integer DataArray, 1 = stage 1, 2 = stage 2") : nothing
    !all(map(x -> x in (1, 2), data[2])) ? error("second column of 'data' can only hold integers 1 or 2 (stages)") : nothing
    n1plan          = interimsamplesize(design)
    n1obs::Integer  = data[1][data[2] .== 1] |> length
    x1obs::Integer  = data[1][data[2] .== 1] |> x -> x[!DataArrays.isna(x)] |> sum
    k::Integer      = data[1][data[2] .== 1][1:min(n1plan, n1obs)] |> x -> x[!DataArrays.isna(x)] |> sum
    nna1::Integer   = data[1][data[2] .== 1][1:min(n1plan, n1obs)] |> x -> x[DataArrays.isna(x)]  |> length
    nobs::Integer   = size(data, 1)
    # adapt to stage one
    stage_one_adapted = design
    if (n1obs != n1plan) | (nna1 != 0) # need to recalculate
        push!(design.params.samplespace.specialnvalues, n1obs) # guarantee that observed n1 is in sample space
        ipm = IPModel(samplespace(parameters(design)), n1obs)
        completemodel(ipm, parameters(design), n1obs)
        for kk in k:(k + nna1)
            addconditionaltypeoneerrorrateconstraint!(ipm, n1obs, kk, nna1, design)
        end
        addinvarianceimputationconstraint!(ipm, n1obs, x1obs, nna1)
        @constraint(ipm.m, getvariable(ipm.m, :xi) <= 0)
        # obj = getobjective(ipm.m)
        # @objective(ipm.m, Min, obj + eta*getvariable(ipm.m, :xi))
        setsolver(ipm.m, solver)
        status = solve(ipm.m)
        if !(status in (:Optimal, :UserLimit)) # no valid solution found!
            warn("no feasible solution reached when adapting to stage one, relaxing conditional type one error rate constraint")
            ipm = IPModel(samplespace(parameters(design)), n1obs)
            completemodel(ipm, parameters(design), n1obs)
            for kk in k:(k + nna1)
                addconditionaltypeoneerrorrateconstraint!(ipm, n1obs, kk, nna1, design)
            end
            addinvarianceimputationconstraint!(ipm, n1obs, x1obs, nna1)
            # @constraint(ipm.m, getvariable(ipm.m, :xi) <= eta)
            obj = getobjective(ipm.m)
            @objective(ipm.m, Min, obj + eta*getvariable(ipm.m, :xi))
            setsolver(ipm.m, solver)
            status = solve(ipm.m)
            if !(status in (:Optimal, :UserLimit)) # no valid solution found!
                error("no feasible solution reached")
            end
            println(getvalue(getvariable(ipm.m, :xi)))
        end
        stage_one_adapted = extractsolution(ipm, parameters(design))
    end
    stage_two_adapted = stage_one_adapted
    n2obs::Integer  = data[1][data[2] .== 2] |> length
    if (n2obs != samplesize(stage_one_adapted, x1obs) - interimsamplesize(stage_one_adapted)) & (n2obs > 0) # only recalculate design after stage two if stage two has positive enrollement
        x2obs::Integer  = data[1][data[2] .== 2] |> x -> x[!DataArrays.isna(x)] |> sum
        l::Integer      = data[1][data[2] .== 2][1:min(n2obs, samplesize(stage_one_adapted, x1obs))] |> x -> x[!DataArrays.isna(x)] |> sum
        nna1            = data[1][data[2] .== 1] |> x -> x[DataArrays.isna(x)]  |> length
        nna2::Integer   = data[1][data[2] .== 2][1:min(n2obs, samplesize(stage_one_adapted, x1obs))] |> x -> x[DataArrays.isna(x)]  |> length
        ipm = IPModel(samplespace(parameters(stage_one_adapted)), n1obs)
        completemodel(ipm, parameters(stage_one_adapted), n1obs)
        # for kk in k:(k + nna1)
        #     addconditionaltypeoneerrorrateconstraint!(ipm, n1obs, kk, nna1, stage_one_adapted)
        # end
        addinvarianceimputationconstraint!(ipm, n1obs, x1obs, nna1)
        addconditionaltypeoneerrorrateconstraint_stagetwo!(ipm, n1obs, x1obs, nna1, l, nna2, stage_one_adapted)
        cvals =  [-Inf; getcvals(samplespace(parameters(stage_one_adapted)), n1obs); Inf]
        for x1 in x1obs:(x1obs + nna1) # for all possible stage 1 outcome n must equal observed value
            @constraint(ipm.m,
                sum(ipm.y[x1, n1obs + n2obs, c] for c in cvals) >= 1
            )
        end
        @constraint(ipm.m, getvariable(ipm.m, :xi) <= 0)
        setsolver(ipm.m, solver)
        status = solve(ipm.m)
        if !(status in (:Optimal, :UserLimit)) # no valid solution found!
            warn("no feasible solution reached when adapting to stage two, relaxing conditional type one error rate control")
            ipm = IPModel(samplespace(parameters(stage_one_adapted)), n1obs)
            completemodel(ipm, parameters(stage_one_adapted), n1obs)
            # for kk in k:(k + nna1)
            #     addconditionaltypeoneerrorrateconstraint!(ipm, n1obs, kk, nna1, stage_one_adapted)
            # end
            addinvarianceimputationconstraint!(ipm, n1obs, x1obs, nna1)
            addconditionaltypeoneerrorrateconstraint_stagetwo!(ipm, n1obs, x1obs, nna1, l, nna2, stage_one_adapted)
            cvals =  [-Inf; getcvals(samplespace(parameters(stage_one_adapted)), n1obs); Inf]
            for x1 in x1obs:(x1obs + nna1) # for all possible stage 1 outcome n must equal observed value
                @constraint(ipm.m,
                    sum(ipm.y[x1, n1obs + n2obs, c] for c in cvals) >= 1
                )
            end
            obj = getobjective(ipm.m)
            @objective(ipm.m, Min, obj + eta*getvariable(ipm.m, :xi))
            setsolver(ipm.m, solver)
            status = solve(ipm.m)
            if !(status in (:Optimal, :UserLimit)) # no valid solution found!
                println(stage_one_adapted)
                println(x1obs)
                println(n2obs)
                println(l)
                error("no feasible solution reached")
            end
            println(getvalue(getvariable(ipm.m, :xi)))
        end
        stage_two_adapted = extractsolution(ipm, parameters(stage_one_adapted))
    end
    # # TODO: check whether constraints also hold for p < p0!
    return stage_two_adapted
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
    @variable(ipm.m, xi >= 0)
    if n1old <= n1obs
        for x1 in x1obs:(x1obs + nna1)
            @constraint(ipm.m,
                sum(
                    dbinom(j, n1obs - n1old, p0)*_cpr(x1 + j, n1obs, n, c, p0)*ipm.y[x1 + j, n, c] for
                    j = 0:(n1obs - n1old),
                    n = nvals,
                    c = cvals
                ) <= _cpr(x1, n1old, samplesize(design, x1), criticalvalue(design, x1), p0) + xi
            )
        end
    end
    if n1old > n1obs
        for x1 in x1obs:(x1obs + nna1)
            @constraint(ipm.m,
                sum(
                    _cpr(x1, n1obs, n, c, null(params))*ipm.y[x1, n, c] for
                    n = nvals,
                    c = cvals
                ) <= sum([dbinom(j, n1old - n1obs, p0)*_cpr(x1 + j, n1old, samplesize(design, x1 + j), criticalvalue(design, x1 + j), p0) for
                    j in 0:(n1old - n1obs)]) + xi
            )
        end
    end
    return ipm
end
function addconditionaltypeoneerrorrateconstraint_stagetwo!(
    ipm, n1, x1obs, nna1, l, nna2, design
)
    params = parameters(design)
    ss     = samplespace(params)
    n1old  = interimsamplesize(design)
    p0     = design |> parameters |> null
    nmax   = maxsamplesize(ss, n1)
    nvals  = getnvals(ss, n1)
    cvals =  [-Inf; getcvals(ss, n1); Inf]
    @variable(ipm.m, xi >= 0)
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
                return 1.0
            end
            if c_x1 == Inf
                return 0.0
            end
            return 1 - Distributions.cdf(Distributions.Binomial(n_x1 - samplesize(design, x1), p0), c_x1 - x1 - l) # TODO: replace
        end
    end
    # rhs
    function rhs(x1, l, n_x1)
        if l > min(n_x1, samplesize(design, x1)) - n1 # not possible
            return 0.0
        end
        if n_x1 < samplesize(design, x1) # new design has smaller value
            # catch c = +- Inf
            if samplesize(design, x1) == -Inf
                return 1.0
            end
            if criticalvalue(design, x1) == Inf
                return 0.0
            end
            return 1 - Distributions.cdf(Distributions.Binomial(samplesize(design, x1) - n_x1, p0), criticalvalue(design, x1) - x1 - l) # TODO: replace
        else
            if x1 + l <= criticalvalue(design, x1)
                return 0.0
            else
                return 1.0
            end
        end
    end
    for x1 in x1obs:(x1obs + nna1)
        for ll in l:(l + nna2)
            @constraint(ipm.m,
                sum(ipm.y[x1, n, c]*(lhs(x1, ll, n, c) - rhs(x1, ll, n)) for
                    n in nvals,
                    c in cvals
                ) <= getvariable(ipm.m, :xi)
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
                    ) - sum(ipm.y[x1 - 1, n, c] for
                        c = cvals
                    ) == 0
                )
            end
        end
    end
    return ipm
end

# function adapt{T<:Integer, TS<:MathProgBase.AbstractMathProgSolver}(
#     design::BinaryTwoStageDesign,
#     outcome::DataArrays.DataVector{Bool},
#     stage::Vector{T},
#     solver::TS
# )
#     length(stage) == length(outcome) ? nothing : throw(InexactError())
#     # check that parameters is not NoParameters (no optimization routine!)
#     params = parameters(design)
#     typeof(params) != NoParameters ? nothing : throw(InexactError())
#     # check that stage is only 1, 2 (cannot be missing)
#     all([s in [1 2] for s in stage]) ? nothing : throw(InexactError())
#     length(stage) == length(outcome) ? nothing : throw(InexactError())
#     # check that parameters is not NoParameters (no optimization routine!)
#     params = parameters(design)
#     typeof(params) != NoParameters ? nothing : throw(InexactError())
#     # check that stage is only 1, 2 (cannot be missing)
#     all([s in [1 2] for s in stage]) ? nothing : throw(InexactError())
#     n1planned = interimsamplesize(design)
#     nmax      = maxsamplesize(params)
#     n1obs     = length(stage[stage .== 1])
#     n1min     = min(n1planned, n1obs)
#     x1obs     = outcome[stage .== 1] |> x -> x[!DataArrays.isna(x)] |> sum
#     nna1      = sum(DataArrays.isna(outcome[stage .== 1][1:n1min]))
#     xobs      = outcome |> x -> x[!DataArrays.isna(x)] |> sum
#     nna       = sum(DataArrays.isna(outcome))
#     nobs      = length(stage)
#     if (n1planned == n1obs) & (nobs == samplesize(design, x1obs)) & (nna == 0) # nothing to do
#         return design
#     end
#     # recreate the same problem conditional on the observed stage one sample size
#     m, y = _createProblem(n1obs, params)
#     m, y = _addconditionaltypeonestageoneconstraint(m, y, n1obs, x1obs, nna1, design)
#     m, y = _addinvarianceimputationstageoneconstraint(m, y, n1obs, x1obs, nna1, design)
#     if any(stage .== 2)
#         m, y = _addconditionaltypeonestagetwoconstraint(m, y, n1planned, x1obs, nna1, xobs, nna, design)
#         cvalsfinite = 0:(nmax - 1)
#         if allowsstoppingforefficacy(params)
#             cvals = [-Inf; cvalsfinite; Inf]
#         else
#             cvals = [cvalsfinite; Inf]
#         end
#         println(x1obs:(x1obs + nna1))
#         for x1 in x1obs:(x1obs + nna1)
#             @constraint(m,
#                 sum(y[x1, nobs, c] for c in cvals) >= 1
#             )
#         end
#     end
#     setsolver(m, solver)
#     status = solve(m)
#     if !(status in (:Optimal, :UserLimit)) # no valid solution found!
#         error("no feasible solution reached")
#     end
#     # extract solution
#     newDesign = _extractSolution(y, n1obs, params) # c.f. util.jl
#     # # TODO: check whether constraints also hold for p < p0!
#     return newDesign
# end
