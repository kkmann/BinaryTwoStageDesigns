type IPModel
    y # the binary assignment variables
    m # the JuMP model
    nvals
    cvals
    function IPModel(ss::SampleSpace, n1::Integer)
        n1 > 0 ? nothing
        ss  = samplespace(params)
        !possible(n1, ss) ? throw(InexactError()) : nothing
        nmax = maxsamplesize(ss, n1)
        nvals = getnvals(ss, n1)
        cvalsfinite = 0:(maximum(nvals) - 1)
        if allowsstoppingforefficacy(params)
            cvals = [-Inf; cvalsfinite; Inf]
            cvalsinfinite = [-Inf; Inf]
        else
            cvals = [cvalsfinite; Inf]
            cvalsinfinite = [Inf]
        end
    end
