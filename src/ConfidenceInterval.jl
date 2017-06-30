abstract ConfidenceInterval

# these two enable array broadcasting of all methods!
Base.size(::ConfidenceInterval) = ()
Base.getindex(ci::ConfidenceInterval, i) = ci

limits{T<:Integer}(ci::ConfidenceInterval, x1::T, x2::T) = error("not implemented!")

confidence(ci::ConfidenceInterval) = try return ci.confidence catch error("not implemented") end

design(ci::ConfidenceInterval) = try return ci.design catch error("not implemented") end


function coverage{T<:Real}(ci::ConfidenceInterval, p::T; orientation = "overall") # or upper or lower
    checkp(p)
    supp   = support(design(ci))
    res    = 0.0
    for i in 1:size(supp, 1)
        x1, x2 = supp[i, :]
        lim = limits(ci, x1, x2)
        if (orientation == "overall") & (lim[1] <= p) & (lim[2] >= p)
            res += pdf(design(ci), x1, x2, p)
        elseif (orientation == "upper") & (lim[2] >= p)
            res += pdf(design(ci), x1, x2, p)
        elseif (orientation == "lower") & (lim[1] <= p)
            res += pdf(design(ci), x1, x2, p)
        end
    end
    return res
end

function meanwidth{T<:Real}(
    ci::ConfidenceInterval,
    p::T
)
    d = design(ci)
    supp   = support(d)
    res    = 0.0
    for i in 1:size(supp, 1)
        x1, x2 = supp[i, :]
        res = res + pdf(d, x1, x2, p)*(limits(ci, x1, x2)[2] - limits(ci, x1, x2)[1])
    end
    return res
end

function meaninterval{T<:Real}(
    ci::ConfidenceInterval,
    p::T
)
    d = design(ci)
    supp   = support(d)
    mean_limits = [0.0; 0.0]
    for i in 1:size(supp, 1)
        x1, x2      = supp[i, :]
        mean_limits = mean_limits + pdf(d, x1, x2, p)*limits(ci, x1, x2)
    end
    return mean_limits
end

function findinconsistencies{T<:Real}(
    ci::ConfidenceInterval,
    p0::T
)
    d = design(ci)
    supp   = support(d)
    res    = []
    for i in 1:1:size(supp, 1)
        x1, x2 = supp[i, :]
        if (limits(ci, x1, x2)[1] <= p0) & (x1 + x2 > d.c[x1 + 1])
            push!(res, [x1 x2 d.c[x1 + 1] limits(ci, x1, x2)[1] p0])
        end
        if (limits(ci, x1, x2)[1] > p0) & (x1 + x2 <= d.c[x1 + 1])
            push!(res, [x1 x2 d.c[x1 + 1] limits(ci, x1, x2)[1] p0])
        end
    end
    res = vcat(res...)
    return res
end
