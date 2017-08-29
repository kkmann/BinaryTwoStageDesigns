"""
    Design{T1<:Integer, T2<:Real, PType<:Parameters}

Binary two-stage design object defining a discrete sampe size function 
``n(x_1)`` and ``c(x_1)``, where ``x_1`` is the observed number of stage-one events 
and the test rejects the superiority null hypothesis whenever ``x_1 + x_2 > c(x_1)``.
"""
mutable struct Design{T1<:Integer, T2<:Real, PType<:Parameters}
    
  n::Vector{T1}
  c::Vector{T2} # rejected if x1 + x2 > c(x1), must allow real to hold Infinity
  params::PType

  function Design{T1,T2,PType}(
    n::Vector{T1}, c::Vector{T2}, params::PType
  ) where {T1<:Integer, T2<:Real, PType<:Parameters}

    any(n .< (length(n) - 1)) ? 
      error("n must always be larger than n1 (= length(n) - 1") : nothing
    length(n) != length(c) ? 
      error("n and c must be of equal length") : nothing
    new(n, c, params)

  end

end # Design


function Design(
  n::Vector{T1}, c::Vector{T2}, params::PType
) where {T1<:Integer,T2<:Real,PType<:Parameters}
  
  return Design{T1,T2,PType}(n, c, params)

end


"""
    Design(n, c, params::PType) where {PType<:Parameters}

Construct a [`Design`](@ref) object from `n` and `c` directly. 
Here, `n` must be convertable to an integer vector and `c` to a real vector 
(both of same length ``n_1``). 
Then, `n[x1 + 1]` is the final sample size after observing `x1` responses in
stage one and `c[x1 + 1]` the corresponding critical value.
A constructor without explicitly specifying a set of parameters with respect to
which thee design is optimal also exists.

# Parameters

| Parameter  | Description |
| ---------: | :---------- |
| n          | sample size function (must be convertable to integer vector) |
| c          | critical value function (must be convertable to real vector) |
| params     | parameter object |

"""
function Design(n, c, params::PType) where {PType<:Parameters}
  
  n = convert(Vector{Int}, n)
  c = convert(Vector{Float64}, c)
  return Design{Int,Float64,Ptype}(n, c, params)

end


"""
    Design(n, c)

Construct a [`Design`](@ref) object from `n` and `c` directly. 
Here, `n` must be convertable to an integer vector and `c` to a real vector 
(both of same length ``n_1``). 

# Parameters

| Parameter  | Description |
| ---------: | :---------- |
| n          | sample size function (must be convertable to integer vector) |
| c          | critical value function (must be convertable to real vector) |

"""
function Design(n, c)
  
  n = convert(Vector{Int}, n)
  c = convert(Vector{Float64}, c)
  return Design{Int,Float64,NoParameters}(n, c, NoParameters())

end


"""
    convert(::Type{DataFrames.DataFrame}, design::Design)

Convert a design to a `DataFrame` with column names `x1`, `n`, and `c`.
"""
convert(::Type{DataFrames.DataFrame}, design::Design) =
    DataFrames.DataFrame(
        x1 = 0:interimsamplesize(design),
        n  = samplesize(design),
        c  = criticalvalue(design)
    )

# these enable array broadcasting of all methods!
Base.size(::Design) = ()

Base.length(design::Design) = 1

Base.start(design::Design) = true

Base.done(design::Design, state) = true

Base.getindex(design::Design, i) = design

Base.show(io::IO, design::Design) = print("Design")


"""
    interimsamplesize(design::Design)

Return interim sample size ``n_1`` of a design.
"""
interimsamplesize(design::Design) = length(design.n) - 1


"""
    parameters(design::Design)

Return stored parameter object of a design.
"""
parameters(design::Design) = design.params


"""
    samplesize(design::Design)

Return vector of final sample sizes.
"""
samplesize(design::Design) = design.n

"""
    samplesize(design::Design, x1::T) where {T<:Integer}

Return final sample size for `x1`.
"""
samplesize(design::Design, x1::T) where {T<:Integer} =
    (checkx1(x1, design); design.n[x1 + 1])


"""
    criticalvalue(design::Design)

Return vector of critical values.
"""
criticalvalue(design::Design) = design.c

"""
    criticalvalue(design::Design, x1::T) where {T<:Integer}

Return critical value for `x1`.
"""
criticalvalue(design::Design, x1::T) where {T<:Integer} =
    (checkx1(x1, design); design.c[x1 + 1])

"""
    pdf{T1<:Integer, T2<:Real}(design::Design, x1::T1, x2::T1, p::T2)

Probability density function of ``(x1, x2)`` responses in stage one and two, 
respectively, under `design` given response rate ``p``.

# Parameters

| Parameter    | Description |
| -----------: | :---------- |
| design       | a Design |
| x1           | number of stage-one responses |
| x2           | number of stage-two responses |
| p            | response probability |

"""
function pdf(
  design::Design,
  x1::T1, x2::T1, p::T2
) where {T1<:Integer, T2<:Real}

  @checkprob p
  try
    checkx1x2(x1, x2, design)
  catch
    return 0.0
  end
  n1 = interimsamplesize(design)
  n  = samplesize(design, x1)
  try
    return p^(x1 + x2)*(1 - p)^(n - x1 - x2)*binomial(n1, x1)*binomial(n - n1, x2) # is 0 if impossible
  catch
    return p^(x1 + x2)*(1 - p)^(n - x1 - x2)*binomial(BigInt(n1), BigInt(x1))*binomial(BigInt(n - n1), BigInt(x2)) # is 0 if impossible
  end

end


"""
    power(design::Design, x1::T1, p::T2) where {T1<:Integer, T2<:Real}

conditional power of a design for a given response rate ``p`` and
the observed number of stage-one responses `x1`.

# Parameters

| Parameter    | Description |
| -----------: | :---------- |
| design       | a [`Design`](@ref) |
| p            | response probability |
| x1           | number of stage-one responses |

"""
function power(
  design::Design, x1::T1, p::T2
) where {T1<:Integer, T2<:Real}
  
  @checkprob p
  checkx1(x1, design)
  res = _cpr(
    x1,
    interimsamplesize(design),
    samplesize(design, x1),
    criticalvalue(design, x1),
    p
  )
  return min(1, max(0, res)) # guarantee bounds!

end


"""
    power(design::Design, p::T) where {T<:Real}

Power of a design for a given response rate ``p`` .

# Parameters

| Parameter    | Description |
| -----------: | :---------- |
| design       | a [`Design`](@ref) |
| p            | response probability |

"""
function power(design::Design, p::T) where {T<:Real}
    
  @checkprob p
  n1      = interimsamplesize(design)
  X1      = Distributions.Binomial(n1, p) # stage one responses
  x1range = collect(0:n1)
  return min(1, max(0, vecdot(Distributions.pdf(X1, x1range), power.(design, x1range, p))))

end


"""
    expectedpower(
      design::Design, x1::T, prior::Function; mcrv::Real = mcrv(parameters(design))
    ) where {T<:Integer}

Compute the conditional expected power of a design given `x1`.

# Parameters

| Parameter    | Description |
| -----------: | :---------- |
| design       | a Design |
| x1           | stage one number of responses |
| prior        | prior function prior(p) for response probability p |
| mcrv         | minimal clinically relevant value for expected power calculation |

"""
function expectedpower(
    design::Design, x1::T, prior::Function; mcrv::Real = mcrv(parameters(design))
) where {T<:Integer}

  checkx1(x1, design)
  z   = QuadGK.quadgk(
      p -> prior(p),  # f(p)
      mcrv,           # p_min
      1,              # p_max
      abstol = 0.001  # tolerance
  )[1]
  res = QuadGK.quadgk(
      p -> prior(p)*power(design, x1, p)/z, # f(p)
      mcrv,         # p_min
      1,             # p_max
      abstol = 0.001 # tolerance
  )[1]
  return min(1, max(0, res)) # guarantee bounds!

end


"""
    expectedpower(design::Design, prior::Function; mcrv::Real = mcrv(parameters(design)))

Compute the expected power of a given design.

# Parameters

| Parameter    | Description |
| -----------: | :---------- |
| design       | a Design |
| prior        | prior function prior(p) for response probability p |

"""
function expectedpower(design::Design, prior::Function; mcrv::Real = mcrv(parameters(design)) )

  z   = QuadGK.quadgk(
      p -> prior(p),             # f(p)
      mcrv,  # p_min
      1,                         # p_max
      abstol = 0.001             # tolerance
  )[1]
  res = QuadGK.quadgk(
      p -> prior(p)*power(design, p)/z, # f(p)
      mcrv,    # p_min
      1,                           # p_max
      abstol = 0.001               # tolerance
  )[1]
  return min(1, max(0, res)) # guarantee bounds!

end



"""
    stoppingforfutility{T<:Real}(design::Design, p::T) where {T<:Real}

Compute probability of stopping early for futility of a given design.

# Parameters

| Parameter    | Description |
| -----------: | :---------- |
| design       | a [`Design`](@ref) |
| p            | response probability |

"""
function stoppingforfutility(design::Design, p::T) where {T<:Real}

  @checkprob p
  n1  = interimsamplesize(design)
  X1  = Distributions.Binomial(n1, p) # stage one responses
  res = 0.0
  c   = criticalvalue(design)
  for x1 in 0:n1
    if c[x1 + 1] == Inf
      res += Distributions.pdf(X1, x1)
    end
  end
  return res

end


"""
    stoppingforefficacy{T<:Real}(design::Design, p::T) where {T<:Real}

Compute probability of stopping early for futility of a given design.

# Parameters

| Parameter    | Description |
| -----------: | :---------- |
| design       | a [`Design`](@ref) |
| p            | response probability |

"""
function stoppingforefficacy(design::Design, p::T) where {T<:Real}

  @checkprob p
  n1  = interimsamplesize(design)
  X1  = Distributions.Binomial(n1, p) # stage one responses
  res = 0.0
  c   = criticalvalue(design)
  for x1 in 0:n1
    if c[x1 + 1] == -Inf
      res += Distributions.pdf(X1, x1)
    end
  end
  return res

end


"""
    test(design::Design, x1::T, x2::T)::Bool where {T<:Integer}

Binary test decision of `design` when observing `x1` responses in stage one and
`x2` responses in stage two.

# Parameters

| Parameter    | Description |
| -----------: | :---------- |
| design       | a [`Design`](@ref) |
| x1           | number of stage-one responses |
| x2           | number of stage-two responses |

"""
function test(design::Design, x1::T, x2::T)::Bool where {T<:Integer}

  checkx1x2(x1, x2, design)
  return x1 + x2 > criticalvalue(design, x1)

end


"""
    simulate(design::Design, p::T2, nsim::T1) where {T1<:Integer, T2<:Real}

Simulate `nsim` trials of design `design` with true response probability `p`.

# Parameters

| Parameter    | Description |
| -----------: | :---------- |
| design       | a [`Design`](@ref) |
| p            | true response probability |
| nsim         | number of simulated trials |

# Return Value

A `DataFrame` with columns `x1` (stage one responses) `n` (overall sample size)
`c` (critical value), `x2` (stage-two responses), and `rejectedH0` (test decision).

"""
function simulate(design::Design, p::T2, nsim::T1) where {T1<:Integer, T2<:Real}

  x2    = SharedArray{Int}(nsim)
  n     = SharedArray{Int}(nsim)
  c     = SharedArray{Float64}(nsim)
  rej   = SharedArray{Bool}(nsim)
  n1    = interimsamplesize(design)
  rv_x1 = Distributions.Binomial(n1, p)
  x1    = SharedArray{Int}(nsim)
  @sync @parallel for i in 1:nsim
      x1[i]  = rand(rv_x1)
      n2     = samplesize(design, x1[i]) - n1
      n[i]   = n2 + n1
      c[i]   = criticalvalue(design, x1[i])
      x2[i]  = rand(Distributions.Binomial(n2, p))
      rej[i] = test(design, x1[i], x2[i])
  end
  return DataFrames.DataFrame(
      x1 = convert(Vector{Int}, x1),
      n  = convert(Vector{Int}, n),
      c  = convert(Vector{Float64}, c),
      x2 = convert(Vector{Int}, x2),
      rejectedH0 = convert(Vector{Bool}, rej)
  )

end


"""
    score(design::Design, params::Parameters)::Real

Evaluates the score of a design for a given set of parameters.

# Parameters

| Parameter    | Description |
| -----------: | :---------- |
| design       | a [`Design`](@ref) |
| params       | a parameters object holding the required information about the score |

# Details

A method `score(design::Design)::Real` is also implemented which
uses the parameter object stored within the design object after optimization.
Note that this is only possible if the design was created as result of calling
`getoptimaldesign()`.

"""
score(design::Design, params::Parameters)::Real = error("not implemented")

score(design::Design)::Real = score(design, parameters(design))


"""
    jeffreysprior(design::Design)

Computes the Jeffreys prior of any given design.

# Parameters

| Parameter    | Description |
| -----------: | :---------- |
| design       | a Design |

# Return Value

A function with signature `prior{T<:Real}(p::T)::Real` where `p` is the response
probability and `prior(p)` the PDF of the Jeffres prior at `p`.

# Examples
```julia-repl
julia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)
julia> params = SimpleMinimalExpectedSampleSize(ss, .2, .4, .05, .2, .4)
julia> design, res = getoptimaldesign(params, solver = Gurobi.GurobiSolver())
julia> f = jeffreysprior(design)
julia> f(.5)
```
"""
function jeffreysprior(design::Design)

  function sqrtfi(p::Float64)::Float64

    res  = 0.0
    n1   = interimsamplesize(design)
    nmax = maximum(samplesize(design))
    supp = support(design)
    for i in 1:size(supp, 1)
        x1, x2 = supp[i, 1], supp[i, 2]
        x   = x1 + x2
        n   = samplesize(design, x1)
        res += binomial(BigInt(n1), BigInt(x1))*binomial(BigInt(n - n1), BigInt(x2))*BigFloat(p^x*(1 - p)^(n - x)*(x/p - (n - x)/(1 - p))^2)
    end
    return sqrt(res)
    
  end
  
  z = QuadGK.quadgk(sqrtfi, 0, 1, abstol = 0.001)[1] # exact integration from 0 to 1 is expensive!
  
  function prior{T<:Real}(p::T)::Real
      
    @checkprob p
    return sqrtfi(p)/z

  end

  return prior

end


"""
    save(filename::String, design::Design)

Save design as .jls file.

# Parameters

| Parameter    | Description |
| -----------: | :---------- |
| filename     | filename for sving design |
| design       | a [`Design`](@ref) |

"""
function save(filename::String, design::Design)
  
  file = open(filename, "w")
  serialize(file, design)
  close(file)

end


"""
    writecsv(filename::String, design::Design; label::String = "") 

Save design as .csv file.

# Parameters

| Parameter    | Description |
| -----------: | :---------- |
| filename     | filename for sving design |
| design       | a [`Design`](@ref) |
| label        | string name of design |
"""
function writecsv(filename::String, design::Design; label::String = "") 

  df = DataFrames.DataFrame(design)
  label != "" ? df[:design] = label : nothing
  DataFrames.writetable(filename, df)
  return df

end


function writepropertiescsv(
  filename::String, design::Design; 
  label::String = "", pvec = linspace(0, 1, 101), nmax::Int = Int(round(1.1*maximum(samplesize(design))))
) 
  
  df = DataFrame(
    p      = Real[],
    n      = Integer[],
    PDF    = Real[],
    ESS    = Real[],
    Q1SS   = Real[],
    Q3SS   = Real[],
    Power  = Real[]
  )
  for p in pvec
    ss = SampleSize(design, p)
    for n in 0:nmax
      push!(df, [p n pdf(ss, n ) mean(ss) quantile(ss, .25) quantile(ss, .75) power(design, p)])
    end
  end
  label != "" ? df[:design] = label : nothing
  DataFrames.writetable(filename, df)
  return df

end




# utility functions
function checkx1(x1::T, design::Design) where {T<:Integer}

  x1 < 0 ? throw(InexactError("x1 must be non-negative")) : nothing
  x1 > interimsamplesize(design) ? throw(InexactError("x1 smaller or equal to n1")) : nothing

end

function checkx1x2(x1::T, x2::T, design::Design) where {T<:Integer}

  checkx1(x1, design)
  x2 < 0 ? throw(InexactError("x2 must be non-negative")) : nothing
  n1 = interimsamplesize(design)
  n2 = samplesize(design, x1) - n1
  x2 > n2 ? throw(InexactError("x2 must be smaller or equal to n2")) : nothing

end

function ispossible(design::Design, x1::T, x2::T) where {T<:Integer}
  
  res = x1 < 0 ? false : true
  res = x1 > interimsamplesize(design) ? false : true
  res = x2 < 0 ? false : true
  res = x2 > samplesize(design, x1) - interimsamplesize(design) ? false : true

end

function support(design::Design)

    n1     = interimsamplesize(design)
    nmax   = maximum(samplesize(design))
    return [[x1, x2] for x1 in 0:n1, x2 in 0:(nmax - n1) if
        (x2 <= samplesize(design, x1) - n1) & ((samplesize(design, x1) > n1) | (x2 == 0))
    ] |> x-> hcat(x...)'

end

function _cpr(x1, n1, n, c, p)
    # conditional probability to reject
    if x1 > c
        return(1.0)
    end
    if n - n1 + x1 <= c
        return(0.0)
    end # TODO: simplify
    return 1 - Distributions.cdf(Distributions.Binomial(n - n1, p), convert(Int64, c - x1))
end

function _cprprior(x1, n1, n, c, p0, prior)
    if x1 > c
        return 1.0
    end
    if n - n1 + x1 <= c
        return 0.0
    end
    z = QuadGK.quadgk(p -> prior(p)*dbinom(c - x1, n - n1, p), p0, 1.0, abstol = 1e-4)[1]
    cposterior(p) = prior(p)*dbinom(c - x1, n - n1, p)/z
    return QuadGK.quadgk(p -> cposterior(p)*_cpr(x1, n1, n, c, p), p0, 1.0, abstol = 1e-4)[1]
end
