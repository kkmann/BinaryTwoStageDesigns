function getoptimaldesign{T<:Integer, TS<:MathProgBase.AbstractMathProgSolver}(
    n1::T,
    parameters::Parameters,
    solver::TS;
    VERBOSE::Integer = 0
)

    !possible(n1, samplespace(parameters)) ? error("n1 not compatible with sample space") : nothing
    # define problem
    # m, y = _createProblem(n1, parameters)
    ipm = IPModel(samplespace(parameters), n1)
    completemodel(ipm, parameters, n1)
    setsolver(ipm.m, solver)
    solve(ipm.m) in (:Optimal, :UserLimit) ? nothing : error("no feasible solution reached")
    design = extractsolution(ipm, parameters)
    _isfeasible(design, parameters) ? nothing : error("solution infeasible")
    return design
end

function getoptimaldesign{TS<:MathProgBase.AbstractMathProgSolver}(
    parameters::Parameters,
    solver::TS;
    VERBOSE::Integer = 1
)
    n1range   = samplespace(parameters) |> interimsamplesizerange
    designs   = Array{BinaryTwoStageDesign}(length(n1range))
    scores    = zeros(length(n1range))
    for i in 1:length(n1range)
        VERBOSE > 0 ? println(n1range[i]) : nothing
        VERBOSE > 0 ? tic() : nothing
        try
            designs[i] = getoptimaldesign(n1range[i], parameters, solver, VERBOSE = VERBOSE)
            scores[i]  = score(designs[i])
        catch e
            println(e)
            scores[i] = Inf
        end
        if VERBOSE > 0
            toc()
            println(scores[i])
            try print(stairs(collect(0:n1range[i])/n1range[i], designs[i].n, title = "n vs. x1/n1", xlim = [0, 1], canvas = AsciiCanvas, style = :pre))
            catch e
                println(e)
            end
            try print(stairs(n1range[1:i], scores[1:i], title = "score vs. n1", xlim = [minimum(n1range), maximum(n1range)], canvas = AsciiCanvas, style = :pre))
            catch e
                println(e)
            end
        end
    end
    return designs[findmin(scores)[2]], Dict(
        "n1"      => n1range,
        "scores"  => scores,
        "designs" => designs
    )
end
