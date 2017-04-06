# function _addconditionaltypeonestageoneconstraint(
#     m, y, n1obs, x1obs, nna1, design::BinaryTwoStageDesign
# )
#     params = parameters(design)
#     ss     = samplespace(params)
#     nvals  = getnvals(ss, n1)
#     if ss.stepsize != 1
#         error("stepsize must be one")
#     end
#     p0     = null(params)
#     n1old  = interimsamplesize(design)
#     nmax   = maxsamplesize(params, n1)
#     cvalsfinite = 0:(maximum(nvals) - 1)
#     if allowsstoppingforefficacy(params)
#         cvals = [-Inf; cvalsfinite; Inf]
#         cvalsinfinite = [-Inf; Inf]
#     else
#         cvals = [cvalsfinite; Inf]
#         cvalsinfinite = [Inf]
#     end
#     if n1old < n1obs
#         for x1 in x1obs:(x1obs + nna1)
#             @constraint(m,
#                 sum(
#                     dbinom(j, n1obs - n1old, p0)*_cpr(x1 + j, n1obs, n, c, p0)*y[x1 + j, n, c] for
#                     j = 0:(n1obs - n1old),
#                     n = n1obs:nmax,
#                     c = cvals
#                 ) <= _cpr(x1, n1old, samplesize(design, x1), criticalvalue(design, x1), p0)
#             )
#         end
#     end
#     if n1old > n1obs
#         for x1 in x1obs:(x1obs + nna1)
#             @constraint(m,
#                 sum(
#                     _cpr(x1, n1obs, n, c, null(params))*y[x1, n, c] for
#                     n = n1obs:nmax,
#                     c = [-Inf; 0:(nmax - 1); Inf]
#                 ) <= sum([dbinom(j, n1old - n1obs, p0)*_cpr(x1 + j, n1old, samplesize(design, x1 + j), criticalvalue(design, x1 + j), p0) for
#                     j in 0:(n1old - n1obs)])
#             )
#         end
#     end
#     return m, y
# end
# function _addconditionaltypeonestagetwoconstraint(
#     m, y, n1, x1obs, nna1, xobs, nna, design
# )
#     params = parameters(design)
#     ss     = samplespace(params)
#     nvals  = getnvals(ss, n1)
#     if ss.stepsize != 1
#         error("stepsize must be one")
#     end
#     p0     = null(params)
#     n1old  = interimsamplesize(design)
#     nmax   = maxsamplesize(params, n1)
#     nvals  = n1old:nmax
#     cvalsfinite = 0:(maximum(nvals) - 1)
#     if allowsstoppingforefficacy(params)
#         cvals = [-Inf; cvalsfinite; Inf]
#         cvalsinfinite = [-Inf; Inf]
#     else
#         cvals = [cvalsfinite; Inf]
#         cvalsinfinite = [Inf]
#     end
#     function lhs(x1::Int64, l::Int64, n_x1::Int64, c_x1::Float64) # c must accept float for +- Inf
#         if l > min(n_x1, samplesize(design, x1)) - n1 # not possible
#             return 0.0
#         end
#         if n_x1 <= samplesize(design, x1) # new design has smaller value
#             if x1 + l <= c_x1
#                 return 0.0
#             else
#                 return 1.0
#             end
#         else
#             # catch c = +- Inf
#             if c_x1 == -Inf
#                 return(1.0)
#             end
#             if c_x1 == Inf
#                 return(0.0)
#             end
#             return 1 - Distributions.cdf(Distributions.Binomial(n_x1 - samplesize(design, x1), p0), c_x1 - x1 - l) # TODO: replace
#         end
#     end
#     # rhs
#     function rhs(x1::Int64, l::Int64, n_x1::Int64, c_x1::Float64)
#         if l > min(n_x1, samplesize(design, x1)) - n1 # not possible
#             return 0.0
#         end
#         if n_x1 < samplesize(design, x1) # new design has smaller value
#             # catch c = +- Inf
#             if samplesize(design, x1) == -Inf
#                 return(1.0)
#             end
#             if criticalvalue(design, x1) == Inf
#                 return(0.0)
#             end
#             return 1 - Distributions.cdf(Distributions.Binomial(samplesize(design, x1) - n_x1, p0), samplesize(design, x1) - x1 - l) # TODO: replace
#         else
#             if x1 + l <= samplesize(design, x1)
#                 return 0.0
#             else
#                 return 1.0
#             end
#         end
#     end
#     for x1 in x1obs:(x1obs + nna1)
#         for l in (xobs - x1obs):(xobs - x1obs + nna - nna1)
#             @constraint(m,
#                 sum(y[x1, n, c]*(lhs(x1, l, n, c) - rhs(x1, l, n, c)) for
#                     n in nvals,
#                     c in cvals
#                 ) <= 0.0
#             )
#         end
#     end
#     return m, y
# end
#
# function _addinvarianceimputationstageoneconstraint(
#     m, y, n1obs, x1obs, nna1, design
# )
#     params = parameters(design)
#     ss     = samplespace(params)
#     nvals  = getnvals(ss, n1)
#     if ss.stepsize != 1
#         error("stepsize must be one")
#     end
#     nmax   = maxsamplesize(params, n1)
#     cvalsfinite = 0:(maximum(nvals) - 1)
#     if allowsstoppingforefficacy(params)
#         cvals = [-Inf; cvalsfinite; Inf]
#         cvalsinfinite = [-Inf; Inf]
#     else
#         cvals = [cvalsfinite; Inf]
#         cvalsinfinite = [Inf]
#     end
#     if nna1 > 0
#         for x1 in (x1obs + 1):(x1obs + nna1)
#             for n in n1obs:nmax
#                 @constraint(m,
#                     sum(y[x1, n, c] for
#                         c = cvals
#                     ) - sum(y[x1 - 1, n, c] for
#                         c = cvals
#                     ) == 0
#                 )
#             end
#         end
#     end
#     return m, y
# end

dbinom(k, n, p) = Distributions.pdf(Distributions.Binomial(n, p), k)

qnorm(p) = Distributions.quantile(Distributions.Normal(0, 1), p)

function findgrid(prior, l, u, n; resolution = 10000)
    candidates = collect(linspace(l, u, resolution + 2))[2:(resolution + 1)]
    quantiles  = collect(linspace(l, u, n + 2))[2:(n + 1)]
    modprior(p) = .9*prior(p) +  .1*Distributions.pdf(Distributions.Beta(.5, .5), p)
    cdfrelaxed = cumsum(modprior.(candidates))
    cdfrelaxed = cdfrelaxed / cdfrelaxed[resolution]
    cdf = cumsum(prior.(candidates))
    cdf = cdf / cdf[resolution]
    pivots = zeros(n)
    cdfpiv = zeros(n)
    i = 1
    for j in 1:n
        while (cdfrelaxed[i] < quantiles[j]) & (i < resolution)
            i += 1
        end
        pivots[j] = candidates[i]
        cdfpiv[j] = cdf[i]
    end
    pivots = (pivots .+ [l; pivots[1:(n - 1)]] ) / 2 # midpoint pivots!
    dcdf = cdfpiv - [0; cdfpiv[1:(n - 1)]] # vector of first order differences
    dcdf = dcdf/sum(dcdf)
    return pivots, dcdf
end

macro checkprob(x...)
    tmp = Expr(:block)
    for xx in x
        name = string(xx)
        push!(tmp.args, :($xx < 0 ? error(@sprintf("%s must be >= 0, is %.5f", $name, $xx)) : true))
        push!(tmp.args, :($xx > 1 ? error(@sprintf("%s must be <= 1, is %.5f", $name, $xx)) : true))
    end
    return(tmp)
end
