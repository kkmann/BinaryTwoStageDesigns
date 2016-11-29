immutable StageOneAdaptedOptimalBinaryTwoStageDesign <: AbstractBinaryTwoStageDesign
    n
    c
    parameters
    originalDesign
    n1obs
    x1obs
    nna1
    function StageOneAdaptedOptimalBinaryTwoStageDesign(oldDesign::AbstractBinaryTwoStageDesign, newDesign::AbstractBinaryTwoStageDesign, n1obs, x1obs, nna1)
        new(newDesign.n, newDesign.c, oldDesign.parameters, oldDesign, n1obs, x1obs, nna1
        )
    end
end

function adaptToInterim(design::OptimalBinaryTwoStageDesign, stageOneData::DataArrays.DataVector{Bool}; CONDITIONALPOWER = true) # isnt x1obs!
    n1old = getInterimSampleSize(design)
    n1obs = length(stageOneData)
    @assert n1obs <= maximum(design.parameters.n1range)
    n1MinOldObs = min(n1old, n1obs)
    x1obs = stageOneData[1:n1MinOldObs] |> x -> x[~DataArrays.isna(x)] |> sum
    nna1  = sum(DataArrays.isna(stageOneData[1:n1MinOldObs]))
    # recreate the same problem conditional on the observed stage one sample size
    m, y, n1obs, params = _createProblem(n1obs, design.parameters)
    # case distinction for n1 <=> n1observed
    if (n1old == n1obs) & (nna1 == 0)
        # nothing to do
        return(design)
    end
    m, y = _addConditionalTypeOneErrorRateStageOneConstraint(m, y, n1obs, x1obs, nna1, design, design.parameters)
    if CONDITIONALPOWER
        m, y = _addConditionalTypeTwoErrorRateStageOneConstraint(m, y, n1obs, x1obs, nna1, design, design.parameters)
    end
    m, y = _addInvarianceUnderImputationStageOneConstraint(m, y, n1obs, x1obs, nna1, design, design.parameters)
    setsolver(m, design.parameters.solver)
    status = solve(m)
    if !(status in (:Optimal, :UserLimit)) # no valid solution found!
        error("no feasible solution reached")
    end
    # extract solution
    newDesign = _extractSolution(y, n1obs, design.parameters.nmax) # c.f. util.jl
    # TODO: check whether constraints also hold for p < p0!
    return StageOneAdaptedOptimalBinaryTwoStageDesign(design, newDesign, n1obs, x1obs, nna1)
end

function adaptToFinal{T<:Integer}(design::BinaryTwoStageDesign, x1obs::T, nna1::T, nobs::T, xobs::T, nna::T, params::PointAlternative)
    # for brevity, extract planning parameters:
    nmax   = params.nmax
    alpha  = params.alpha
    beta   = params.beta
    p0     = params.p0
    p1     = params.p1
    solver = params.solver
    @assert nobs <= nmax
    n1 = getInterimSampleSize(design)
    # recreate initial problem
    m, y, n1, params = _createProblem(n1, params)
    # case distinction for n <=> nobs
    if nobs == getSampleSize(design)
        # nothing to do
        return(design)
    end
    # na constraints: x1obs .. x1obs + nna1 must have equal n!
    if nna1 > 0
        for x1 in (x1obs + 1):(x1obs + nna1)
            for n in n1:nmax
                @constraint(m,
                    sum{y[x1, n, c], c = [-Inf; 0:(nmax - 1); Inf]} - sum{y[x1 - 1, n, c], c = [-Inf; 0:(nmax - 1); Inf]} == 0
                )
            end
        end
    end
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
            return 1 - Distributions.cdf(Distributions.Binomial(n_x1 - design.n[x1 + 1], p0), c_x1 - x1 - l) # TODO: replace
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
            return 1 - Distributions.cdf(Distributions.Binomial(design.n[x1 + 1] - n_x1, p0), design.c[x1 + 1] - x1 - l) # TODO: replace
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
    # add constraints that ensures n(x1obs) = nobs (and for all other possible outcomes of stage one)
    for x1 in x1obs:(x1obs + nna1)
        @constraint(m,
            sum{y[x1, nobs, c], c = [-Inf; 0:(nmax -1); Inf]} >= 1
        )
    end
    setsolver(m, solver)
    status = solve(m)
    if !(status in (:Optimal, :UserLimit)) # no valid solution found!
        error("no feasible solution reached")
    end
    # get solution
    design = _extractSolution(y, n1, params.nmax)
    return design
end
