type IPModel
    y # the binary assignment variables
    m # the JuMP model
    nvals
    cvalsfinite
    cvals
    n1
    ss
    function IPModel(ss::SampleSpace, n1::Integer)
        n1 <= 1           ? error("n1 must be >= 1")                     : nothing
        !possible(n1, ss) ? error("n1 not compatible with sample space") : nothing
        nmax = maxsamplesize(ss, n1)
        nvals = getnvals(ss, n1)
        cvalsfinite = getcvals(ss, n1)
        cvalsinfinite = [-Inf; Inf]
        cvals =  [-Inf; cvalsfinite; Inf]
        m = Model()
        # indicator variables y[x1, n, c] == 1 iff n(x1) = n, c(x1) = c
        @variable(m, y[x1 = 0:n1, n = nvals, c = cvals], Bin)
        if !isgroupsequential(ss) # dummy variables for unimodality
            @variable(m, _is_mode[x1 = 0:n1], Bin)
        else
            @variable(m, cont[x1 = 0:n1], Bin) # dummy variables for determining
                                               # whether design continues
        end
        for x1 in 0:n1
            if isgroupsequential(ss)
                @constraint(m, # lhs is one if design continues at x1
                    sum(y[x1, n, c] for n in nvals, c in cvalsinfinite) - (1 - cont[x1]) == 0
                )
            end
            @constraint(m, # functional constraint: exactly one non-zero entry in (y) per x1
                sum(y[x1, n, c] for n in nvals, c in cvals) == 1
            )
            @constraint(m, # functional constraint: exactly one non-zero entry in (y) per x1
                n*sum(y[x1, n, c] for n in nvals, c in cvals) >= n1
            )
            if x1 > 0 # contiguous stopping for futility
                @constraint(m,
                    y[x1 - 1, n1, Inf] - y[x1, n1, Inf] >= 0
                )
            end
            if (x1 < n1)
                @constraint(m, # contiguous stopping for efficacy
                    y[x1 + 1, n1, -Inf] - y[x1, n1, -Inf] >= 0
                )
            end
            for n in nvals
                if isgroupsequential(ss)
                    @constraint(m, # groupsequential designs can only have fixed decision with early stopping
                        sum(y[x1, n1, c] for c in cvalsinfinite) == sum(y[x1, n, c] for n in nvals, c in cvalsinfinite)
                    )
                end
                for c in cvals
                    if isfinite(c)
                        if c >= n
                            @constraint(m, y[x1, n, c] == 0)
                        end
                    end
                    if isfinite(c) & (n == n1) # c finite => n > n1
                        @constraint(m, y[x1, n, c] == 0)
                    end
                    if (x1 > c) & (c != -Inf) # decision is fixed but c is not valid
                        @constraint(m, y[x1, n, c] == 0)
                    end
                    if !possible(n1, n, c, ss) # sample space constraints
                        @constraint(m,
                            y[x1, n, c] == 0
                        )
                    end
                    if isgroupsequential(ss) & isfinite(c)
                        if x1 < n1
                            @constraint(m, # n, c must be equal for all x1 with finite c
                                y[x1, n, c] - y[x1 + 1, n, c] - (2 - cont[x1] - cont[x1 + 1]) <= 0
                            )
                            @constraint(m, # n, c must be equal for all x1 with finite c
                                y[x1, n, c] - y[x1 + 1, n, c] + (2 - cont[x1] - cont[x1 + 1]) >= 0
                            )
                        end
                    end
                end
            end
            if !isgroupsequential(ss)
                for x1_ in 1:x1
                    @constraint(m,
                        sum(n*(y[x1_, n, c] - y[x1_ - 1, n, c]) for
                            n in nvals[2:end],
                            c in cvals
                        ) - 3*nmax*_is_mode[x1] >= -3*nmax
                    )
                end
                for x1_ in (x1 + 1):n1
                    @constraint(m,
                        sum(n*(y[x1_, n, c] - y[x1_ - 1, n, c]) for
                            n in nvals[2:end],
                            c in cvals
                        ) + 3*nmax*_is_mode[x1] <= 3*nmax
                    )
                end
            end
        end
        if !isgroupsequential(ss)
            @constraint(m, # unimodality
                sum(_is_mode[x1] for x1 in 0:n1) >= 1 # at least one mode (can be on boundary as well!)
            )
        end
        new(y, m, nvals, cvalsfinite, cvals, n1, ss)
    end
end

getnvals(ipm::IPModel, n1::Integer) = getnvals(ipm.ss, n1)
getcvals(ipm::IPModel, n1::Integer) = getcvals(ipm.ss, n1)

function extractsolution(ipm::IPModel, params)
    nvec = zeros(Int64, ipm.n1 + 1)
    cvec = zeros(Float64, ipm.n1 + 1) # need float for +/- Inf
    val::Real = 0.0
    current_best_score::Real = 0.0
    for x1 in 0:ipm.n1
        current_best_score = 0.0
        for n in ipm.nvals
            for c in ipm.cvals
                val = getvalue(ipm.y[x1, n, c])
                if val > current_best_score # in the end gives the solution where y[x1, n, c] is closest to 1
                    current_best_score = val
                    nvec[x1 + 1] = n
                    cvec[x1 + 1] = c
                end
            end
        end
        if current_best_score < .5 # too far from integrality, set to zero and try other fix
            nvec[x1 + 1] = 0
            cvec[x1 + 1] = 0
        end
    end
    try
        return BinaryTwoStageDesign(nvec, cvec, params)
    catch e
        println(nvec)
        println(cvec)
        try
            for i in 1:length(nvec)
                if nvec[i] == 0
                    if (i > 1) & !isfinite(cvec[i - 1])
                        nvec[i] = nvec[i - 1]
                        cvec[i] = cvec[i - 1]
                        continue
                    end
                    if (i < length(nvec)) & !isfinite(cvec[i + 1])
                        nvec[i] = nvec[i + 1]
                        cvec[i] = cvec[i + 1]
                    end
                end
            end
            println(nvec)
            println(cvec)
            design = BinaryTwoStageDesign(nvec, cvec, params)
            return design
        catch e
            rethrow(e)
        end
    end
end
