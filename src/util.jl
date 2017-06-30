dbinom(k, n, p) = Distributions.pdf(Distributions.Binomial(n, p), k)

qnorm(p) = Distributions.quantile(Distributions.Normal(0, 1), p)

function findgrid(prior, l, u, n)
    dp = 1/(n + 1)
    modprior_(p) = .9*prior(p) + .1*Distributions.pdf(Distributions.Beta(1, 1), p)
    zmodprior = QuadGK.quadgk(modprior_, l, u)[1]
    modprior(p) = modprior_(p) / zmodprior * (p < l ? 0 : 1)
    zprior = QuadGK.quadgk(prior, l, u)[1]
    modcdf(p) = QuadGK.quadgk(modprior, l, p)[1]
    cprior(p) = prior(p) / zprior * (p < l ? 0 : 1)
    cdf(p) = QuadGK.quad(cprior, l, p)[1]
    pivots  = [l + 10^-6.0]
    modcdfs = [0.0]
    cdfs    = [0.0]
    for i in 2:(n + 1)
        push!(pivots, Roots.fzero(p -> modcdf(p) - dp - modcdfs[i - 1], pivots[i - 1], 1))
        push!(modcdfs, modcdf(pivots[i]))
        push!(cdfs, cdf(pivots[i]))
    end
    pivots = (pivots[2:(n + 1)] .+ pivots[1:n] ) / 2 # midpoint pivots!
    dcdf = cdfs[2:(n + 1)] - cdfs[1:n] # vector of first order differences
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
