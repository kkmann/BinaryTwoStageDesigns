"""
    optimaldesign(
      n1::Integer,
      parameters::Parameters,
      solver::MathProgBase.AbstractMathProgSolver
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

Design object optimized for given parameter set.

# Examples
```julia-repl
julia> ss     = SimpleSampleSpace(10:25, 100, n2min = 5)
julia> design = optimaldesign(15, params, solver = Gurobi.GurobiSolver())
```
"""
function optimaldesign(
    n1::Integer,
    parameters::Parameters,
    solver::MathProgBase.AbstractMathProgSolver
)

    !possible(n1, samplespace(parameters)) ? warn("n1 not compatible with sample space") : nothing
    ipm = IPModel(samplespace(parameters), n1, UNIMODAL = isunimodal(parameters))
    completemodel(ipm, parameters, n1)
    JuMP.setsolver(ipm.m, solver)
    JuMP.solve(ipm.m) in (:Optimal, :UserLimit) ? nothing : error("no feasible solution reached")
    design = extractsolution(ipm, parameters)
    _isfeasible(design, parameters) ? nothing : error("solution infeasible")
    return design

end # optimaldesign


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

Tuple, first element is the Design object optimized for given parameter set
and the second element is a dictionary with n1, scores and respective optimal designs.

# Examples
```julia-repl
julia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)
julia> interimsamplesize(ss)
julia> design, res = getoptimaldesign(params, solver = Gurobi.GurobiSolver())
```
"""
function optimaldesign(
  parameters::Parameters,
  solver::MathProgBase.AbstractMathProgSolver;
  VERBOSE::Integer = 1,
  EARLYTERMINATION::Bool = false
)

  n1range   = samplespace(parameters) |> interimsamplesizerange
  designs   = Array{Design}(length(n1range))
  scores    = ones(length(n1range))*Inf

  if VERBOSE > 0
    t0 = now()
    println()

    println(@sprintf(
      "optimizing design for parameters '%s'", label(parameters)
    ))
    print(@sprintf(
      "\rconsidering %i stage-one sample sizes between %i and %i", length(n1range), minimum(n1range), maximum(n1range)
    ))
    print(@sprintf(
      " using %s as solver\n\r", typeof(solver)
    ))
  end

  function reldiff(newval, best)  
    
    diff = best - newval
    if diff == 0
      return 0.0
    else 
      return -sign(diff) * abs(diff/best)
    end

  end

  for i in 1:length(n1range)
    t0i = now()
    try
      designs[i] = optimaldesign(n1range[i], parameters, solver)
      scores[i]  = score(designs[i])
      if VERBOSE > 0
        println()
        println("    time    n1   % done   sol. time [s]   cum. time [min]       score        best   % diff to best")
        print(@sprintf("\r%s   %3i    %5.1f   %13i   %15.1f   %+.2e   %+.2e   %14.1f", 
          Dates.format(t0i, "HH:MM:SS"), 
          n1range[i], 
          i / length(n1range) * 100, 
          Int(round(Int(Dates.value(now() - t0i)) / 1000)), 
          Int(Dates.value(now() - t0)) / 60 / 1000, 
          scores[i],
          scores[findmin(scores)[2]],
          reldiff(scores[i], scores[findmin(scores)[2]]) * 100
        ))
        println()
      end
    catch e
      println()
      print(@sprintf("\r%s   %3i   error: %s, continuing...", Dates.format(t0i, "HH:MM:SS"), n1range[i], e))
      println()
      scores[i] = Inf
    end
    if VERBOSE > 1
      try 
        println()
        print(UnicodePlots.stairs(
          collect(0:n1range[i]) / n1range[i], 
          designs[i].n, 
          title = "n vs. x1/n1", 
          ylim = [0, Int(round(maximum(samplesize(designs[i])*1.1)))],
          canvas = UnicodePlots.AsciiCanvas, 
          style = :pre
        ))
        print("\n\n")
      catch e
      end
      println(convert(DataFrame(design[i])))
    end
    print("\n")
    if (EARLYTERMINATION == true) & (i > 3)
      try
        if all(reldiff.(scores[[i - 2; i - 1; i]], scores[findmin(scores)[2]]) .>= .05)
          info("terminating optimization early")
          break
        end
      catch e
      end
    end
  end
  print("\n\n")
  if VERBOSE > 1
    valid = scores .!= Inf
    try 
      print(UnicodePlots.stairs(n1range[valid], scores[valid], title = "score vs. n1", xlim = [minimum(n1range), maximum(n1range)], canvas = UnicodePlots.AsciiCanvas, style = :pre))
      print("\n\n")
    catch e
    end
  end
  if VERBOSE > 0
    print(@sprintf("\rdone after %i minutes.\n\r", Int(round(Int(Dates.value(now() - t0)) / 60 / 1000))))
    print(@sprintf("\rOptimal stage-one sample size is %i resulting in a minimal score of %.2e\n\n\r", n1range[findmin(scores)[2]], scores[findmin(scores)[2]]))
  end
  return designs[findmin(scores)[2]], Dict("n1" => n1range, "scores" => scores, "designs" => designs)

end
