immutable OptimalBinaryTwoStageDesign <: AbstractBinaryTwoStageDesign
    n
    c
    parameters
    function OptimalBinaryTwoStageDesign(parameters::Parameters)
        design, res = optimalDesign(parameters)
        new(
            design.n,
            design.c,
            parameters
        )
    end
    function OptimalBinaryTwoStageDesign(design, parameters::Parameters)
        new(
            design.n,
            design.c,
            parameters
        )
    end
end

function optimalDesign(
    n1::Int64,      # stage one sample size
    parameters::Parameters
)
    @assert n1 in parameters.n1range
    # define problem
    m, y, n1, parameters = _createProblem(n1, parameters)
    setsolver(m, parameters.solver)
    # solve
    status = solve(m)
    if !(status in (:Optimal, :UserLimit)) # no valid solution found!
        error("no feasible solution reached")
    end
    # extract solution
    design = _extractSolution(y, n1, parameters.nmax) # c.f. util.jl
    return OptimalBinaryTwoStageDesign(design, parameters)
end

function optimalDesign(parameters::Parameters)
    designs = pmap(
        n1 -> try
            optimalDesign(n1, parameters)
        catch e
            println(e)
        end,
        parameters.n1range
    )
    scores = map(design -> try parameters.score(design) catch e Inf end, designs) # if score cannot be computed, probably infeasible > Inf
    return designs[findmin(scores)[2]], Dict(
        "n1"      => collect(parameters.n1range),
        "scores"  => scores,
        "designs" => designs
    )
end
