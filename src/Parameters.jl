abstract Parameters # only guarantees p0, nmax, n1range, alpha # TODO: also score, expectedScore


# make parameters iterable for automatic broadcasting
Base.size(::Parameters) = ()
Base.getindex(par::Parameters, i) = par


null(par::Parameters) = try par.p0 catch error("not implemented") end
alpha(par::Parameters) = try par.alpha catch error("not implemented") end


samplespace(par::Parameters) = try par.samplespace catch error("not implemeted") end
maxsamplesize(par::Parameters) = error("not implemented")
regularization(par::Parameters) = error("not implemented")
efficacy(par::Parameters) = error("not implemented")




type NoParameters <: Parameters
end
null(par::NoParameters) = error("NoParameters do not have null hypothesis")
alpha(par::NoParameters) = error("NoParameters do not have alpha")
samplespace(par::NoParameters) = error("NoParameters do not have sample space")
maxsamplesize(par::NoParameters) = error("NoParameters do not have maximal sample size")




abstract PointAlternative <: Parameters
alternative(par::PointAlternative) = try par.p1 catch error("not implemented") end
beta(par::PointAlternative) = try par.beta catch error("not implemented") end




abstract VagueAlternative <: Parameters
prior{T<:Real}(par::VagueAlternative, p::T) = error("not implemented")
