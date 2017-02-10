type CompatibleEstimator <: BinaryTwoStageDesignEstimator
    design::AbstractBinaryTwoStageDesign
    estimates::DataFrames.DataFrame

    function CompatibleEstimator(                 # for any practial case it will
        design::AbstractBinaryTwoStageDesign,     # be much more efficient to
        solver;                                   # precompute all estimates in
        prior::Function = jeffreysprior(design),  # advance!
        k = 100
    )
        epsilon = 10.0^(-6.0)
        isodd(k) ? throw(InexactError()) : nothing
        n1   = interimsamplesize(design)
        nmax = maximum(samplesize(design))
        supp = support(design)
        model = Model(solver = solver)
        @variable(model, 0.0 <= est[x1 = 0:n1, x2 = 0:(nmax - n1); ispossible(design, x1, x2)] <= 1.0)
        h = 1/k
        vrho = collect(linspace(0, 1, k + 1))
        function coefs(x1, x2, rho)
            if !ispossible(design, x1, x2)
                return 0.0
            else
                pdf = prior(rho)
                pdf = (pdf > 10.0) | !isfinite(pdf) ? 10 : pdf
                return pdf*rho^(x1 + x2)*(1.0 - rho)^(samplesize(design, x1) - x1 - x2)
            end
        end
        @objective(model, Min,
            sum(binomial(BigInt(n1), BigInt(x1))*binomial(BigInt(samplesize(design, x1) - n1), BigInt(x2))*h/3*(
                      coefs(x1, x2, vrho[2*i - 1])*(est[x1, x2] - vrho[2*i - 1])^2 +
                    4*coefs(x1, x2,     vrho[2*i])*(est[x1, x2] -     vrho[2*i])^2 +
                      coefs(x1, x2, vrho[2*i + 1])*(est[x1, x2] - vrho[2*i + 1])^2) for
                x1 in 0:n1, x2 in 0:(nmax - n1), i in 1:convert(Int64, floor(k/2)) if
                ispossible(design, x1, x2)
            )
        )
        # solve(model)
        # 3.    enforce test consistency by requiring the largest non-rejecting
        #       estimate to be smaller than the smallest rejecting one
        @variable(model, 0 <= minREJECT <= 1)
        @variable(model, 0 <= maxACCEPT <= 1)
        for i in 1:size(supp, 1)
            x1, x2 = supp[i, :]
            if test(design, x1, x2) == true
                @constraint(model, minREJECT - est[x1, x2] <= 0)
            else
                @constraint(model, maxACCEPT - est[x1, x2] >= 0)
            end
        end
        @constraint(model, minREJECT - maxACCEPT - epsilon >= 0) # smallest estimate of reject must be larger than largest estimate of accept
        for i in size(supp, 1)
            for j in size(supp, 1)
                if i < j
                    x1a, x2a = supp[i, :]
                    x1b, x2b = supp[j, :]
                    if (x1a == x1b) & (x2a >= x2b) # same first stage, larger or equal second stage ~~> larger or equal estimate
                        @constraint(model, est[x1a, x2a] - est[x1b, x2b] >= 0)
                    end
                    n2a = samplesize(design, x1a)
                    n2a = samplesize(design, x1a)
                    if (n2a == interimsamplesize(design)) & (n2b == interimsamplesize(design)) & (x1a >= x1b) # no second stage, larger or equal first stage ~~> larger or equal estimate
                        @constraint(model, est[x1a, x2a] - est[x1b, x2b] >= 0)
                    end
                end
            end
        end
        solve(model)
        estimates = Array(Float64, n1 + 1, nmax - n1 + 1)
        for x1 in 0:n1
            for x2 in 0:(nmax - n1) # brute force, not all possible
                if ispossible(design, x1, x2)
                    estimates[x1 + 1, x2 + 1] = getvalue(est[x1, x2])
                else
                    estimates[x1 + 1, x2 + 1] = convert(Float64, NaN)
                end
            end
        end
        new(design, estimates)
    end

end

function estimate{T<:Integer}(estimator::CompatibleEstimator, x1::T, x2::T)
    ispossible(design(estimator), x1, x2) ? nothing : throw(InexactError())
    return estimator.estimates[x1 + 1, x2 + 1]
end
