function getoptimaldesign{T<:Integer, TS<:MathProgBase.AbstractMathProgSolver}(
    n1::T,
    parameters::Parameters,
    solver::TS;
    VERBOSE::Integer = 0
)
    VERBOSE > 0 ? tic() : nothing
    possible(n1, samplespace(parameters)) ? nothing : throw(InexactError())
    # define problem
    m, y = _createProblem(n1, parameters)
    setsolver(m, solver)
    status = solve(m)
    if !(status in (:Optimal, :UserLimit)) # no valid solution found!
        error("no feasible solution reached")
    end
    design = _extractSolution(y, n1, parameters) # c.f. util.jl
    _isfeasible(design, parameters) ? nothing : error("solution infeasible")
    if VERBOSE > 0
        toc()
        println(@sprintf("n1: %i", n1)) 
    end
    return design
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
    return designs[findmin(scores)[2]], Dict(
        "n1"      => n1range,
        "scores"  => scores,
        "designs" => designs
    )
end
