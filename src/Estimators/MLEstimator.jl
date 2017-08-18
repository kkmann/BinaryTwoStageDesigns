"""
    MLEstimator

    MLEstimator(design::Design)

Simple maximum likelihood estimator for response rate `p`.

# Examples
```julia-repl
julia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)
julia> interimsamplesize(ss)
julia> design = getoptimaldesign(15, params, solver = Gurobi.GurobiSolver())
julia> est = MLEstimator(design)
```
"""
struct MLEstimator{TD<:Design} <: Estimator
    
  design::TD

  MLEstimator{TD}(design::TD) where {TD<:Design} = new(design)

end # MLEStimator


MLEstimator(design::TD) where {TD<:Design} = MLEstimator{TD}(design)


Base.show(io::IO, estimator::MLEstimator) = print("MLEstimator")


function estimate(estimator::MLEstimator, x1::T, x2::T) where {T<:Integer}

    checkx1x2(x1, x2, design(estimator))
    return (x1 + x2)/samplesize(design(estimator), x1)

end
