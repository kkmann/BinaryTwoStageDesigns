"""
    Parameters

Abstract type representing a generic set of paramters for finding optimal
two-stage designs. 
"""
abstract type Parameters end


# make parameters iterable for automatic broadcasting
Base.size(::Parameters) = ()

Base.length(par::Parameters) = 1

Base.getindex(par::Parameters, i) = par


"""
    label(par::Parameters)

Return the label of a [`Parameters`](@ref)-object.
"""
label(par::Parameters) = try par.label catch error("not implemented") end

"""
    null(par::Parameters)

Return the response rate under the null hypothesis of a [`Parameters`](@ref)-object.
"""
null(par::Parameters) = try par.p0 catch error("not implemented") end

"""
    mtoer(par::Parameters)

Return the maximal type one error rate of a [`Parameters`](@ref)-object.
"""
mtoer(par::Parameters) = try par.mtoer catch error("not implemented") end

"""
    mcrv(par::Parameters)

Return the minimal clinically relevant response rate of a [`Parameters`](@ref)-object.
"""
mcrv(par::Parameters) = try par.pmcrv catch error("not implemented") end

"""
    samplespace(par::Parameters)

Return the [`SampleSpace`](@ref)-object of a [`Parameters`](@ref)-object.
"""
samplespace(par::Parameters) = try par.samplespace catch error("not implemeted") end

"""
    maxsamplesize(par::Parameters)

Return the maximal overall sample size of a [`Parameters`](@ref)-object.
"""
maxsamplesize(par::Parameters) = error("not implemented")

"""
    isgroupsequential(par::Parameters)::Bool

Returns `true` if the undrlying [`SampleSpace`](@ref) is constraint to 
group-sequential designs.
"""
isgroupsequential(par::Parameters)::Bool = error("not implemented")

"""
    allowsstoppingforfutility(par::Parameters)::Bool

Returns `true` if the undrlying [`SampleSpace`](@ref) allows stopping-for-futlity.
"""
allowsstoppingforfutility(par::Parameters) = error("not implemented")


"""
    NoParameters <: Parameters

This dummy object is used is a [`Design`](@ref) is created manually without 
providing information on the respective optimality criterion.
"""
struct NoParameters <: Parameters end

null(par::NoParameters) = error("NoParameters do not have null hypothesis")

mtoer(par::NoParameters) = error("NoParameters do not have mtoer")

samplespace(par::NoParameters) = error("NoParameters do not have sample space")

maxsamplesize(par::NoParameters) = error("NoParameters do not have maximal sample size")

Base.show(io::IO, object::NoParameters) = print("NoParameters")


"""
    PointAlternative <: Parameters

Abstract type representing a generic set of paramters for finding optimal
two-stage designs which are characterized by a point alternative.
"""
abstract type PointAlternative <: Parameters end

alternative(par::PointAlternative) = try par.p1 catch error("not implemented") end

mtter(par::PointAlternative) = try par.beta catch error("not implemented") end


"""
    VagueAlternative <: Parameters

Abstract type representing a generic set of paramters for finding optimal
two-stage designs which are characterized by a prior distribution over the 
response probability.
"""
abstract type VagueAlternative <: Parameters end

prior(par::VagueAlternative, p::T) where {T<:Real} = try par.prior(p) catch error("not implemented") end
