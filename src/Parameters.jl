"""
    Parameters

Abstract type representing a generic set of paramters for finding optimal
two-stage designs.
"""
abstract Parameters


# make parameters iterable for automatic broadcasting
Base.size(::Parameters) = ()
Base.length(par::Parameters) = 1
Base.getindex(par::Parameters, i) = par


null(par::Parameters) = try par.p0 catch error("not implemented") end
alpha(par::Parameters) = try par.alpha catch error("not implemented") end
mcrv(par::Parameters) = try par.pmcrv catch error("not implemented") end

samplespace(par::Parameters) = try par.samplespace catch error("not implemeted") end
maxsamplesize(par::Parameters) = error("not implemented")
isgroupsequential(par::Parameters) = error("not implemented")
allowsstoppingforfutility(par::Parameters) = error("not implemented")

function simulate(params::Parameters, p, n1, n; nna1 = 0, nna2 = 0, x1 = -1, x2 = -1) # still needed?
    RV = Distributions.Bernoulli(p)
    if x1 < 0
        response1 = convert(DataArrays.DataArray{Bool}, rand(RV, n1))
        response1[Distributions.sample(1:n1, nna1, replace = false)] = DataFrames.NA
    else
        response1 = convert(DataArrays.DataArray{Bool}, fill(false, n1))
        response1[Distributions.sample(1:n1, x1, replace = false)] = true
        negativeresponses = view(response1, response1 .== false)
        negativeresponses[Distributions.sample(1:(n1 - x1), nna1, replace = false)] = DataFrames.NA
    end
    if x2 < 0
        response2 = convert(DataArrays.DataArray{Bool}, rand(RV, n - n1))
        response2[Distributions.sample(1:(n - n1), nna2, replace = false)] = DataFrames.NA
    else
        response2 = convert(DataArrays.DataArray{Bool}, fill(false, n - n1))
        response2[Distributions.sample(1:(n - n1), x2, replace = false)] = true
        negativeresponses = view(response2, response2 .== false)
        negativeresponses[Distributions.sample(1:(n - n1 - x2), nna2, replace = false)] = DataFrames.NA
    end
    data = DataFrames.DataFrame(
        response = [response1; response2],
        stage    = [fill(1, n1); fill(2, n - n1)]
    )
    return data
end

"""
    NoParameters <: Parameters

Empty parameters.
"""
type NoParameters <: Parameters
end
null(par::NoParameters) = error("NoParameters do not have null hypothesis")
alpha(par::NoParameters) = error("NoParameters do not have alpha")
samplespace(par::NoParameters) = error("NoParameters do not have sample space")
maxsamplesize(par::NoParameters) = error("NoParameters do not have maximal sample size")
function show(io::IO, object::NoParameters)
    return "no parameters"
end


"""
    PointAlternative <: Parameters

Abstract type representing a generic set of paramters for finding optimal
two-stage designs which are characterized by a specific point alternative.
"""
abstract PointAlternative <: Parameters
alternative(par::PointAlternative) = try par.p1 catch error("not implemented") end
beta(par::PointAlternative) = try par.beta catch error("not implemented") end


"""
    VagueAlternative <: Parameters

Abstract type representing a generic set of paramters for finding optimal
two-stage designs which are characterized by a prior distribution over the response probability.
"""
abstract VagueAlternative <: Parameters
prior{T<:Real}(par::VagueAlternative, p::T) = try par.prior(p) catch error("not implemented") end
