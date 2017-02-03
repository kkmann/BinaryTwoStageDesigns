abstract Parameters # only guarantees p0, nmax, n1range, alpha # TODO: also score, expectedScore


# make parameters iterable for automatic broadcasting
Base.size(::Parameters) = ()
Base.getindex(par::Parameters, i) = par


function get_null(par::Parameters)::Float64
    try
        return par.p0
    catch
        error("not implemented")
    end
end
function get_n_stage_one(par::Parameters)
    try
        return par.n1range
    catch
        error("not implemented")
    end
end
function get_n_max(par::Parameters)
    try
        return par.nmax
    catch
        error("not implemented")
    end
end
function get_max_type_one_error(par::Parameters)
    try
        return par.alpha
    catch
        error("not implemented")
    end
end

abstract PointAlternative <: Parameters # must also have p1 + beta
function get_alternative(par::PointAlternative)
    try
        return par.p1
    catch
        error("not implemented")
    end
end
function get_min_power(par::PointAlternative)
    try
        return par.beta
    catch
        error("not implemented")
    end
end

abstract VagueAlternative <: Parameters # must implement p0 + Beta prior on p1 but no beta
function get_alternative(par::VagueAlternative)
    try
        return par.prior
    catch
        error("not implemented")
    end
end
function get_min_power(par::VagueAlternative)
    try
        return par.beta
    catch
        error("not implemented")
    end
end
