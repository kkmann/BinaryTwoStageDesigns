function jeffreysPrior(design::AbstractBinaryTwoStageDesign)
    _sqrtFisherInformation = function(p::Float64)
        res  = 0.0
        n1   = getInterimSampleSize(design)
        nmax = maximum(getSampleSize(design))
        supp = _support(design)
        for x1 in 0:n1
            for x2 in 0:(nmax - n1)
                if _isInSupport(supp, x1, x2)
                    x = x1 + x2
                    n = getSampleSize(design, x1)
                    res += binomial(BigInt(n1), BigInt(x1))*binomial(BigInt(n - n1), BigInt(x2))*BigFloat(p^x*(1 - p)^(n - x)*(x/p - (n - x)/(1 - p))^2)
                end
            end
        end
        return sqrt(res)
    end
    constant = quadgk(_sqrtFisherInformation, 0, 1, abstol = 0.01)[1] # exact integration from 0 to 1 is expensive!
    function prior{T<:Real}(p::T)::Float64
        @assert 0 <= p
        @assert p <= 1
        return _sqrtFisherInformation(p)/constant
    end
    return prior
end

type CompatibleEstimator <: BinaryTwoStageDesignEstimator
    design::AbstractBinaryTwoStageDesign
    estimates::DataFrames.DataFrame

    function CompatibleEstimator(design::AbstractBinaryTwoStageDesign, solver; prior::Function = jeffreysPrior(design), k = 100)
        epsilon = 10.0^(-6.0)
        n1   = getInterimSampleSize(design)
        nmax = maximum(getSampleSize(design))
        supp = _support(design)
        ISPOSSIBLE = falses(n1 + 1, nmax - n1 + 1)
        for x1 in 0:n1
            for x2 in 0:(nmax - n1) # brute force, not all possible
                ISPOSSIBLE[x1 + 1, x2 + 1] = _isInSupport(supp, x1, x2)
            end
        end
        possibleInd   = find(x -> x == true, ISPOSSIBLE)
        possible_x1x2 = hcat(ind2sub((n1 + 1, nmax - n1 + 1), possibleInd)...) - 1
        nPossible     = size(possible_x1x2, 1)
        model = Model(solver = solver)
        @variable(model, 0.0 <= est[x1 = 0:n1, x2 = 0:(nmax - n1); ISPOSSIBLE[x1 + 1, x2 + 1]] <= 1.0)
        if isodd(k)
            error("k is not even!")
        end
        h::Float64 = 1/k
        vrho = convert(Vector{Float64}, (0:k)/k)
        # we need to precompute all coefficients as @objective does not accept
        # exponents larger than two (even for not-model variables)
        coefs = Array(Float64, n1 + 1, nmax - n1 + 1, k + 1)
        for i in 1:(k + 1)
            pdf = prior(vrho[i])
            if (pdf > 10.0) | !isfinite(pdf)
                pdf = 10.0
            end
            for x1 in 0:n1
                for x2 in 0:(nmax - n1)
                    if ISPOSSIBLE[x1 + 1, x2 + 1]
                        coefs[x1 + 1, x2 + 1, i] = pdf*vrho[i]^(x1 + x2)*(1.0 - vrho[i])^(design.n[x1 + 1] - x1 - x2)
                    else
                        coefs[x1 + 1, x2 + 1, i] = 0.0
                    end
                end
            end
        end
        @objective(
            model,
            Min,
            sum{ # integrate using composite Simpson's rule!
                binomial(BigInt(n1), BigInt(x1))*binomial(BigInt(convert(Int64, design.n[x1 + 1]) - n1), BigInt(x2))*h/3*(coefs[x1 + 1, x2 + 1, 2*i - 1]*(est[x1, x2] - vrho[2*i - 1])^2 + 4*coefs[x1 + 1, x2 + 1, 2*i]*(est[x1, x2] - vrho[2*i])^2 + coefs[x1 + 1, x2 + 1, 2*i + 1]*(est[x1, x2] - vrho[2*i + 1])^2),
                x1 = 0:n1,
                x2 = 0:(nmax - n1),
                i = 1:convert(Int64, floor(k/2));
                ISPOSSIBLE[x1 + 1, x2 + 1]
            }
        )
        # solve(model)
        # 3.    enforce test consistency by requiring the largest non-rejecting
        #       estimate to be smaller than the smallest rejecting one
        @variable(
            model,
            0 <= minREJECT <= 1
        )
        @variable(
            model,
            0 <= maxACCEPT <= 1
        )
        for x in 1:nPossible
            x1 = possible_x1x2[x, 1]
            x2 = possible_x1x2[x, 2]
            if test(design, x1, x2) == true
                @constraint(
                    model,
                    minREJECT - est[x1, x2] <= 0
                )
            else
                @constraint(
                    model,
                    maxACCEPT - est[x1, x2] >= 0
                )
            end
        end
        @constraint(
            model,
            minREJECT - maxACCEPT - epsilon >= 0 # smallest estimate of reject must be larger than largest estimate of accept
        )

        # 4.    equal first stage
        for t1 in 1:nPossible
            for t2 in 1:nPossible
                if t1 < t2
                    x1a = possible_x1x2[t1, 1]
                    x2a = possible_x1x2[t1, 2]
                    x1b = possible_x1x2[t2, 1]
                    x2b = possible_x1x2[t2, 2]
                    if (x1a == x1b) & (x2a >= x2b)
                        # same first stage, larger or equal second stage ~~> larger or equal estimate
                        @constraint(
                            model,
                            est[x1a, x2a] - est[x1b, x2b] >= 0
                        )
                    end
                end
            end
        end

        # 5.    only first stage
        for t1 in 1:nPossible
            for t2 in 1:nPossible
                if t1 < t2
                    x1a = possible_x1x2[t1, 1]
                    n2a = design.n[x1a + 1]
                    x1b = possible_x1x2[t2, 1]
                    n2b = design.n[x1b + 1]
                    if (n2a == 0) & (n2b == 0) & (x1a >= x1b)
                        # no second stage, larger or equal first stage ~~> larger or equal estimate
                        @constraint(
                            model,
                            est[x1a, x2a] - est[x1b, x2b] >= 0
                        )
                    end
                end
            end
        end
        # 6.    solve the model
        solve(model)

        # 7.    extract solution
        estimates = Array(Float64, n1 + 1, nmax - n1 + 1)
        for x1 in 0:n1
            for x2 in 0:(nmax - n1) # brute force, not all possible
                if ISPOSSIBLE[x1 + 1, x2 + 1]
                    estimates[x1 + 1, x2 + 1] = getvalue(est[x1, x2])
                else
                    estimates[x1 + 1, x2 + 1] = convert(Float64, NaN)
                end
            end
        end

        # 8.    build estimator object
        new(
            design,
            estimates
        )
    end

end

function estimate{T<:Integer}(estimator::CompatibleEstimator, x1::T, x2::T)::Float64
    supp = _support(estimator.design)
    if !_isInSupport(supp, x1, x2)
        return NaN
    else
        return estimator.estimates[x1 + 1, x2 + 1]
    end
end
