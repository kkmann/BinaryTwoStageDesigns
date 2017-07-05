"""
    getoptimaldesign{T<:Integer, TS<:MathProgBase.AbstractMathProgSolver}(
        n1::T,
        parameters::Parameters,
        solver::TS;
        VERBOSE::Integer = 0
    )

Find the optimal two-stage design for given parameters and fixed stage-one
sample size.

# Parameters

| Parameter    | Description |
| -----------: | :---------- |
| n1           | stage-one sample size |
| nparameters  | paramters object defining the optimality criterion |
| solver       | MathProgBase solver used for optimization |
| VERBOSE      | control verbosity during optimization |

# Return Value

BinaryTwoStageDesign object optimized for given parameter set.

# Examples
```julia-repl
julia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)
julia> interimsamplesize(ss)
julia> design = getoptimaldesign(15, params, solver = Gurobi.GurobiSolver())
```
"""
function getoptimaldesign{T<:Integer, TS<:MathProgBase.AbstractMathProgSolver}(
    n1::T,
    parameters::Parameters,
    solver::TS;
    VERBOSE::Integer = 0
)
    !possible(n1, samplespace(parameters)) ? warn("n1 not compatible with sample space") : nothing
    ipm = IPModel(samplespace(parameters), n1)
    completemodel(ipm, parameters, n1)
    JuMP.setsolver(ipm.m, solver)
    JuMP.solve(ipm.m) in (:Optimal, :UserLimit) ? nothing : error("no feasible solution reached")
    design = extractsolution(ipm, parameters)
    _isfeasible(design, parameters) ? nothing : error("solution infeasible")
    return design
end

"""
    getoptimaldesign{TS<:MathProgBase.AbstractMathProgSolver}(
        parameters::Parameters,
        solver::TS;
        VERBOSE::Integer = 1
    )

Find the optimal two-stage design for given parameters (optimizes over `n1` as well).

# Parameters

| Parameter    | Description |
| -----------: | :---------- |
| nparameters  | paramters object defining the optimality criterion |
| solver       | MathProgBase solver used for optimization |
| VERBOSE      | control verbosity during optimization |

# Return Value

Tuple, first element is the BinaryTwoStageDesign object optimized for given parameter set
and the second element is a dictionary with n1, scores and respective optimal designs.

# Examples
```julia-repl
julia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)
julia> interimsamplesize(ss)
julia> design, res = getoptimaldesign(params, solver = Gurobi.GurobiSolver())
```
"""
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
            try print(UnicodePlots.stairs(collect(0:n1range[i])/n1range[i], designs[i].n, title = "n vs. x1/n1", xlim = [0, 1], canvas = UnicodePlots.AsciiCanvas, style = :pre))
            catch e
                println(e)
            end
            try print(UnicodePlots.stairs(n1range[1:i], scores[1:i], title = "score vs. n1", xlim = [minimum(n1range), maximum(n1range)], canvas = UnicodePlots.AsciiCanvas, style = :pre))
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
