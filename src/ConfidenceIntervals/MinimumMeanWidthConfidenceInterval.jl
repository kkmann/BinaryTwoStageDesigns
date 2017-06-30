immutable MinimumMeanWidthConfidenceInterval <: ConfidenceInterval
    estimator::BinaryTwoStageDesignEstimator
    confidence::Float64
    limits::Array{Float64}
    prior::Function

    function MinimumMeanWidthConfidenceInterval{TS<:MathProgBase.AbstractMathProgSolver}(
        estimator::BinaryTwoStageDesignEstimator,
        rho0::Float64,
        prior::Function,
        solver::TS;
        confidence::Float64 = .9,
        ngrid::Int64  = 100
    )
        alpha =  1 - confidence
        ngrid += 1 # number of gridpoints is number of intervals + 1
        d     = design(estimator)
        n1    = interimsamplesize(d)
        nmax::Int64  = maximum(samplesize(d))

        # construct array of possible outcomes
        x1x2_poss = zeros(Int64, 1, 2)
        for x1 in 1:n1
            for x2 in 0:convert(Int64, d.n[x1 + 1] - n1)
                x1x2_poss = [x1x2_poss; convert(Array{Int64, 2}, [x1 x2])]
            end
        end
        nposs = size(x1x2_poss, 1)
        function ispossible(x1, x2)
            sum(mapslices(x -> x == vec([x1 x2]), x1x2_poss, 2)) >= 1
        end

        # set up arrays needed for coefficients
        rho_grid = linspace(0, 1, ngrid)
        # rho_grid_refined = linspace(0, 1, (ngrid - 1)*refinement_fctr)
        probs = zeros(Float64, nposs, ngrid)
        exp_probs = zeros(Float64, nposs)
        estimates = zeros(Int64, nposs, 2)
        for o in 1:nposs
            x1 = x1x2_poss[o, 1]
            x2 = x1x2_poss[o, 2]
            # bracket estimator in grid
            i  = 1
            while rho_grid[i] < estimate(estimator, x1, x2)
                i += 1
            end
            if rho_grid[i] == estimate(estimator, x1, x2)
                estimates[o, 1] = i
                estimates[o, 2] = i
            else
                estimates[o, 1] = i
                estimates[o, 2] = i + 1
            end
            # precompute probabilities
            exp_probs[o] = QuadGK.quadgk(rho -> prior(rho)*pdf(d, x1, x2, rho), 0.0, 1.0, abstol = 10^(-3.0))[1]
            for i in 1:ngrid
                probs[o, i] = pdf(d, x1, x2, rho_grid[i])
            end
        end
        model = JuMP.Model(solver = solver)
        JuMP.@variable(model, # lower bounds
            y[o = 1:nposs, i = 1:(ngrid - 1), j = (i + 1):ngrid],
            Bin
        )
        JuMP.@variable(model,
            z_mon[o = 1:nposs, i = 1:2] >= 0
        )
        a = 0
        for i in 1:ngrid
            if rho0 > rho_grid[i]
                a = i
            end
        end
        for o in 1:nposs
            x1 = x1x2_poss[o, 1]
            x2 = x1x2_poss[o, 2]
            # monotonicity violations
            if o > 1
                x1_prev = x1x2_poss[o - 1, 1]
                x2_prev = x1x2_poss[o - 1, 2]
                if (x1 == x1_prev) & (x2 == x2_prev + 1)
                    expr1 = JuMP.AffExpr()
                    expr2 = JuMP.AffExpr()
                    for i in 1:(ngrid - 1)
                        for j in (i + 1):ngrid
                            JuMP.append!(expr1, rho_grid[i]*(y[o - 1, i, j] - y[o, i, j]))
                            JuMP.append!(expr2, rho_grid[j]*(y[o - 1, i, j] - y[o, i, j]))
                        end
                    end
                    JuMP.@constraint(model,
                        # expr1 <= 0
                        z_mon[o, 1] >= expr1
                    )
                    JuMP.@constraint(model,
                        # expr2 <= 0
                        z_mon[o, 2] >= expr2
                    )
                end
            end
            # functional constraint
            funcConstrExpr = JuMP.AffExpr()
            for i in 1:(ngrid - 1)
                for j in (i + 1):ngrid
                    JuMP.append!(funcConstrExpr, y[o, i, j])
                end
            end
            JuMP.@constraint(model,
                funcConstrExpr == 1
            )
            # estimator consistency: estimator must lie in interval!
            estConstrExpr = JuMP.AffExpr()
            for i in 1:estimates[o, 1] # lower bound must be smaller or equal to estimate
                for j in max(i + 1, estimates[o, 2]):ngrid # upper bound must be larger or equal
                    JuMP.append!(estConstrExpr, y[o, i, j]) # a possible assignment
                end
            end
            JuMP.@constraint(model,
                estConstrExpr == 1 # any of the possible assigments must be selected!
            )
            # consistency 1: x1 + x2 > c(x1) -> lower bound must be strictly larger than rho0
            if x1 + x2 > d.c[x1 + 1]
                cons1ConstrExpr = JuMP.AffExpr()
                for i in 1:min(ngrid - 1, a) # rho_grid[a] <= rho0
                    for j in (i + 1):ngrid
                        JuMP.append!(cons1ConstrExpr, y[o, i, j])
                    end
                end
                JuMP.@constraint(model,
                    cons1ConstrExpr == 0 # none of the assignments with lower bound <= rho0 can be selected
                )
            end
            # consistency 2: x1 + x2 <= c(x1) -> lower bound must be smaller or equal to rho0
            if x1 + x2 <= d.c[x1 + 1]
                cons2ConstrExpr = JuMP.AffExpr()
                for i in 1:min(ngrid - 1, a) # rho_grid[a] <= rho0
                    for j in (i + 1):ngrid
                        JuMP.append!(cons2ConstrExpr, y[o, i, j])
                    end
                end
                JuMP.@constraint(model,
                    cons2ConstrExpr == 1 # one of the assignments with lower bound <= rho0 must be selected
                )
            end
            # upper boundary bleeding, if x1 + x2 == n(x1) the upper boundary must be 1 to prevent coverage from dropping towards zero
            if x1 + x2 == d.n[x1 + 1]
                boundaryUpperConstrExpr = JuMP.AffExpr()
                for i in 1:(ngrid - 1)
                    JuMP.append!(boundaryUpperConstrExpr, y[o, i, ngrid])
                end
                JuMP.@constraint(model,
                    boundaryUpperConstrExpr == 1
                )
            end
        end
        # lower boundary bleeding, if x1 + x2 == 0 lower boundary must be 0 to prevent coverage fro mdropping too low
        boundaryLowerConstrExpr = JuMP.AffExpr()
        for j in 2:ngrid
            JuMP.append!(boundaryLowerConstrExpr, y[1, 1, j]) # x1x2_poss[1, :] = [0 0]
        end
        JuMP.@constraint(model,
            boundaryLowerConstrExpr == 1
        )
        for i in 1:ngrid
            lowerBoundCoverage = JuMP.AffExpr()
            upperBoundCoverage = JuMP.AffExpr()
            for o in 1:nposs
                for j in 1:(ngrid - 1)
                    for k in (j + 1):ngrid
                        if j <= i # lower bound <= rho_grid[i]
                            JuMP.append!(lowerBoundCoverage, probs[o, i]*y[o, j, k])
                        end
                        if k >= i # upper bound >= rho_grid[i]
                            JuMP.append!(upperBoundCoverage, probs[o, i]*y[o, j, k])
                        end
                    end
                end
            end
            JuMP.@constraint(model,
                lowerBoundCoverage >= 1 - alpha/2
            )
            JuMP.@constraint(model,
                upperBoundCoverage >= 1 - alpha/2
            )
        end
        objExpr = JuMP.AffExpr()
        for o in 1:nposs
            JuMP.append!(objExpr, (ngrid - 1)*z_mon[o, 1])
            JuMP.append!(objExpr, (ngrid - 1)*z_mon[o, 2])
            for i in 1:(ngrid - 1)
                for j in (i + 1):ngrid
                    JuMP.append!(objExpr, exp_probs[o]*(rho_grid[j] - rho_grid[i])*y[o, i, j])
                end
            end
        end
        JuMP.@objective(model, Min, objExpr)
        status = JuMP.solve(model)
        if !(status in (:Optimal, :UserLimit)) # no valid solution found!
            error("no feasible solution reached")
        end
        # create actual confidence interval object
        limits = fill(NaN, n1 + 1, nmax - n1 + 1, 2)
        for o in 1:nposs
            x1 = x1x2_poss[o, 1]
            x2 = x1x2_poss[o, 2]
            current_best_score = 0.0
            for i in 1:(ngrid - 1)
                for j in (i + 1):ngrid
                    val = JuMP.getvalue(y[o, i, j]) == 1.0
                    if val > current_best_score # in the end gives the solution where y[x1, n, c] is closest to 1
                        current_best_score = val
                        limits[x1 + 1, x2 + 1, :] = [rho_grid[i] rho_grid[j]]
                    end
                end
            end
        end
        return new(estimator, confidence, limits, prior)
    end
end

estimator(ci::MinimumMeanWidthConfidenceInterval) = ci.estimator

design(ci::MinimumMeanWidthConfidenceInterval) = ci |> estimator |> design


function limits{T<:Integer}(ci::MinimumMeanWidthConfidenceInterval, x1::T, x2::T)
    ispossible(design(ci), x1, x2) ? nothing : throw(InexactError())
    return convert(Vector{Float64}, view([ci.limits[x1 + 1, x2 + 1, 1] ci.limits[x1 + 1, x2 + 1, 2]], 1:2))
end
