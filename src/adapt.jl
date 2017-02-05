function adapt{T<:Integer, TS<:MathProgBase.AbstractMathProgSolver}(
    design::BinaryTwoStageDesign,
    outcome::DataArrays.DataVector{Bool},
    stage::Vector{T},
    solver::TS;
    CONDITIONALPOWER::Bool = true
)
    length(stage) == length(outcome) ? nothing : throw(InexactError())
    # check that parameters is not NoParameters (no optimization routine!)
    params = parameters(design)
    typeof(params) != NoParameters ? nothing : throw(InexactError())
    # check that stage is only 1, 2 (cannot be missing)
    all([s in [1 2] for s in stage]) ? nothing : throw(InexactError())
    n1planned = interimsamplesize(design)
    n1obs     = length(stage[stage .== 1])

    n1min  = min(n1planned, n1obs)
    x1obs  = outcome[stage .== 1] |> x -> x[!DataArrays.isna(x)] |> sum
    nna1   = sum(DataArrays.isna(outcome[stage .== 1][1:n1min]))
    if (n1planned == n1obs) & (nna1 == 0) # nothing to do
        return design
    end
    # recreate the same problem conditional on the observed stage one sample size
    m, y = _createProblem(n1obs, params)
    m, y = _addconditionaltypeonestageoneconstraint(m, y, n1obs, x1obs, nna1, design)
    m, y = _addinvarianceimputationstageoneconstraint(m, y, n1obs, x1obs, nna1, design)
    # if CONDITIONALPOWER
    #     # maintianing conditional type one and type two error rates and
    #     # invariance under imputation is not guaranteed to yield a solution!
    #     # the conditional power constraint is the only optional one which could
    #     # be dropped...
    #     m, y = _addConditionalTypeTwoErrorRateStageOneConstraint(m, y, n1obs, x1obs, nna1, design, design.parameters)
    # end
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
#
# function _adaptstageone(
#     design::StageOneAdaptedOptimalBinaryTwoStageDesign,
#     data::DataFrames.DataVector{Bool},
#     CONDITIONALPOWER = true
# )
#     n1old = getInterimSampleSize(design)
#     n1obs = length(stageOneData)
#     @assert n1obs <= maximum(design.parameters.n1range)
#     n1MinOldObs = min(n1old, n1obs)
#     x1obs = stageOneData[1:n1MinOldObs] |> x -> x[~DataArrays.isna(x)] |> sum
#     nna1  = sum(DataArrays.isna(stageOneData[1:n1MinOldObs]))
#     # recreate the same problem conditional on the observed stage one sample size
#     m, y, n1obs, params = _createProblem(n1obs, design.parameters)
#     if (n1old == n1obs) & (nna1 == 0)
#         # nothing to do
#         return StageOneAdaptedOptimalBinaryTwoStageDesign(design, design, n1obs, x1obs, nna1)
#     end
#     m, y = _addConditionalTypeOneErrorRateStageOneConstraint(m, y, n1obs, x1obs, nna1, design, design.parameters)
#     if CONDITIONALPOWER
#         # maintianing conditional type one and type two error rates and
#         # invariance under imputation is not guaranteed to yield a solution!
#         # the conditional power constraint is the only optional one which could
#         # be dropped...
#         m, y = _addConditionalTypeTwoErrorRateStageOneConstraint(m, y, n1obs, x1obs, nna1, design, design.parameters)
#     end
#     m, y = _addInvarianceUnderImputationStageOneConstraint(m, y, n1obs, x1obs, nna1, design, design.parameters)
#     setsolver(m, design.parameters.solver)
#     status = solve(m)
#     if !(status in (:Optimal, :UserLimit)) # no valid solution found!
#         error("no feasible solution reached")
#     end
#     # extract solution
#     newDesign = _extractSolution(y, n1obs, design.parameters.nmax) # c.f. util.jl
#     # TODO: check whether constraints also hold for p < p0!
#     return newDesign
# end
#
# function adapt(design::StageOneAdaptedOptimalBinaryTwoStageDesign, data::DataFrames.DataFrame)
#     nmax   = design.parameters.nmax
#     solver = design.parameters.solver
#     nobs  = DataFrames.nrow(data)
#     @assert nobs <= nmax
#     n1 = getInterimSampleSize(design)
#     @assert n1 == data |> x -> x[x[:stage] .== 1, :response] |> length
#     @assert :stage in names(data)
#     @assert all(map(x -> x in [1; 2], data[:stage]))
#     @assert :response in names(data)
#     @assert isa(data[:response], DataArrays.DataArray{Bool, 1})
#     x1obs = data |> x -> x[x[:stage] .== 1, :response] |> x -> x[~DataArrays.isna(x)] |> sum
#     nna1  = data |> x -> x[x[:stage] .== 1, :response] |> DataArrays.isna |> sum
#     xobs  = data[:response] |> x -> x[~DataArrays.isna(x)] |> sum
#     nna   = data[:response] |> DataArrays.isna |> sum
#     # recreate initial problem
#     m, y, n1, params = _createProblem(n1, design.parameters)
#     # case distinction for n <=> nobs
#     if (nobs == getSampleSize(design, x1obs)) & (nna == 0)
#         # nothing to do
#         return StageTwoAdaptedOptimalBinaryTwoStageDesign(design, design, data)
#     end
#     m, y = _addConditionalTypeOneErrorRateStageTwoConstraint(m, y, n1, x1obs, nna1, xobs, nna, design, design.parameters)
#     m, y = _addInvarianceUnderImputationStageOneConstraint(m, y, n1, x1obs, nna1, design, design.parameters)
#     # add constraints that ensures n(x1obs) = nobs (and for all other possible outcomes of stage one)
#     for x1 in x1obs:(x1obs + nna1)
#         @constraint(m,
#             sum{y[x1, nobs, c], c = [-Inf; 0:(nmax -1); Inf]} >= 1
#         )
#     end
#     setsolver(m, solver)
#     status = solve(m)
#     if !(status in (:Optimal, :UserLimit)) # no valid solution found!
#         error("no feasible solution reached")
#     end
#     # get solution
#     newDesign = _extractSolution(y, n1, design.parameters.nmax) # c.f. util.jl
#     # TODO: check whether constraints also hold for p < p0!
#     return StageTwoAdaptedOptimalBinaryTwoStageDesign(design, newDesign, data)
# end
