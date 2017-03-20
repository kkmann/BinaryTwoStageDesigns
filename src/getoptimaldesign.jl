function getoptimaldesign{T<:Integer, TS<:MathProgBase.AbstractMathProgSolver}(
    n1::T,
    parameters::Parameters,
    solver::TS;
    VERBOSE::Integer = 0
)
    VERBOSE > 0 ? println(n1) : nothing
    VERBOSE > 0 ? tic() : nothing
    possible(n1, samplespace(parameters)) ? nothing : throw(InexactError())
    # define problem
    m, y = _createProblem(n1, parameters)
    setsolver(m, solver)
    status = solve(m)
    if !(status in (:Optimal, :UserLimit)) # no valid solution found!
        error("no feasible solution reached")
    end
    try
        design = _extractSolution(y, n1, parameters) # c.f. util.jl
        _isfeasible(design, parameters) ? nothing : error("solution infeasible")
        VERBOSE > 0 ? toc() : nothing
        VERBOSE > 0 ? println(status) : nothing
        VERBOSE > 0 ? println(score(design)) : nothing
        VERBOSE > 0 ? println(convert(DataFrames.DataFrame, design)) : nothing
        VERBOSE > 0 ? println() : nothing
        return design
    catch e
         println(e)
         error("could not extract solution")
    end
end

function getoptimaldesign{TS<:MathProgBase.AbstractMathProgSolver}(
    parameters::Parameters,
    solver::TS;
    VERBOSE::Integer = 1
)
    n1range = interimsamplesizerange(samplespace(parameters))
    designs = pmap(
        n1 -> try
            getoptimaldesign(n1, parameters, solver, VERBOSE = VERBOSE)
        catch e
            println(e)
        end,
        n1range
    )
    scores = map(design -> try score(design) catch e Inf end, designs) # if score cannot be computed, probably infeasible > Inf
    VERBOSE > 0 ? println(scores) : nothing
    return designs[findmin(scores)[2]], Dict(
        "n1"      => n1range,
        "scores"  => scores,
        "designs" => designs
    )
end
