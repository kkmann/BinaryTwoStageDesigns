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
