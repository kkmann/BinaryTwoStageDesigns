function getoptimaldesign{T<:Integer, TS<:MathProgBase.AbstractMathProgSolver}(
    n1::T,
    parameters::Parameters,
    solver::TS
)
    possible(n1, samplespace(parameters)) ? nothing : throw(InexactError())
    # define problem
    m, y = _createProblem(n1, parameters)
    setsolver(m, solver)
    status = solve(m)
    if !(status in (:Optimal, :UserLimit)) # no valid solution found!
        error("no feasible solution reached")
    end
    design = _extractSolution(y, n1, parameters) # c.f. util.jl
    return design
end

# function getoptimaldesign(parameters::Parameters)
#     designs = pmap(
#         n1 -> try
#             getoptimaldesign(n1, parameters)
#         catch e
#             println(e)
#         end,
#         interimsamplesizerange(samplespace(parameters))
#     )
#     scores = map(design -> try parameters.score(design) catch e Inf end, designs) # if score cannot be computed, probably infeasible > Inf
#     return designs[findmin(scores)[2]], Dict(
#         "n1"      => collect(parameters.n1range),
#         "scores"  => scores,
#         "designs" => designs
#     )
# end
