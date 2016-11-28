function optimalDesign(
    n1::Int64,      # stage one sample size
    params::PlanningAssumptions
)
    # for brevity, extract planning parameters:
    @assert n1 in params.n1range
    # define problem
    m, y, n1, params = _createProblem(n1, params)
    setsolver(m, params.solver)
    # solve
    status = solve(m)
    if !(status in (:Optimal, :UserLimit)) # no valid solution found!
        error("no feasible solution reached")
    end
    # extract solution
    design = _extractSolution(y, n1, params.nmax) # c.f. util.jl
    return design
end

function optimalDesign(params::PlanningAssumptions)
    designs = pmap(
        n1 -> try
            optimalDesign(n1, params)
        catch e
            println(e)
        end,
        params.n1range
    )
    scores = map(design -> try params.score(design) catch e Inf end, designs) # if score cannot be computed, probably infeasible > Inf
    return designs[findmin(scores)[2]], Dict(
        "n1"      => collect(params.n1range),
        "scores"  => scores,
        "designs" => designs
    )
end

immutable OptimalBinaryTwoStageDesign <: BinaryTwoStageDesign
    n
    c
    parameters
    function BinaryTwoStageDesign(parameters::PlanningAssumptions)
        design, res = optimalDesign(parameters)
        new(design.n, design.c, parameters)
    end
end
