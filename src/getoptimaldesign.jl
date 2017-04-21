function getoptimaldesign{T<:Integer, TS<:MathProgBase.AbstractMathProgSolver}(
    n1::T,
    parameters::Parameters,
    solver::TS;
    VERBOSE::Integer = 0
)
    !possible(n1, samplespace(parameters)) ? warn("n1 not compatible with sample space") : nothing
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
    parameters::Parameters, n1lo, n1hi, solver::TS;
    VERBOSE::Integer = 1
)
    piv::Vector{Int64}   = zeros(0)
    val::Vector{Float64} = zeros(0)
    designs              = Any[]
    function f_(n1::Int)
        ind = findin(piv, n1)
        if length(ind) == 0 # not yet computed
            sval = NaN
            try
                dsg = getoptimaldesign(n1, parameters, solver, VERBOSE = VERBOSE)
                if VERBOSE > 0
                    try print(stairs(collect(0:n1)/n1, dsg.n, title = "n vs. x1/n1", xlim = [0, 1], canvas = AsciiCanvas, style = :pre))
                    catch e
                        println(e)
                    end
                end
                sval = score(dsg)
                println([n1 sval])
                append!(val, sval)
                append!(piv, n1)
                push!(designs, dsg)
            catch e
                println(e)
                sval = Inf
                append!(val, sval)
                append!(piv, n1)
                push!(designs, nothing)
            end
            tmp     = sortperm(piv)
            piv     = piv[tmp]
            val     = val[tmp]
            designs = designs[tmp]
            if VERBOSE > 0
                try print(lineplot(piv[val .!= Inf], val[val .!= Inf], title = "score vs. n1", xlim = [n1lo, n1hi], canvas = AsciiCanvas))
                catch e
                    println(e)
                end
            end
            return sval
        end
        if length(ind) == 1
            return val[ind[1]]
        end
        if !(length(ind) in [0 1])
            error("length(ind) not in [0, 1]")
        end
    end
    f_(n1lo)
    println(val)
    !isfinite(val[1]) ? error("problem must be solvable for n1lo") : nothing
    f_(n1hi)
    !isfinite(val[2]) ? error("problem must be solvable for n1hi") : nothing
    clo, cup = n1lo, n1hi
    while true
        if cup - clo == 2 # only one integer between
            f_(clo + 1)
            break # converged
        end
        if cup - clo == 1 # no integer between, converged
            break
        end
        prop = convert(Vector{Int}, round(linspace(clo, cup, 4))) # proposed pivots, prop[1] and prop[4] already examined
        println(prop)
        # compute new values (or look them up)
        propvals = f_.(prop)
        minind   = findmin(propvals)[2]
        # continue with subinterval containing the minimum
        if minind <= 2 # continue with left subinterval
            clo, cup = prop[1], prop[3]
        else # continue with right
            clo, cup = prop[2], prop[4]
        end
    end
    return designs[findmin(val)[2]], Dict(
        "n1"      => piv,
        "scores"  => val,
        "designs" => designs
    )
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
