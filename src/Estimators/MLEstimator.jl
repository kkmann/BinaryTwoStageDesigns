"""
    MLEstimator

Simple maximum likelihood estimator for response rate `p`.
"""
struct MLEstimator{TD<:Design} <: Estimator
    
  design::TD

  MLEstimator{TD}(design::TD) where {TD<:Design} = new(design)

end # MLEStimator

"""
    MLEstimator(design::TD) where {TD<:Design}

Create maximum likelihood estimator for response rate under given design.
"""
MLEstimator(design::TD) where {TD<:Design} = MLEstimator{TD}(design)


Base.show(io::IO, estimator::MLEstimator) = print("MLEstimator")


function estimate(estimator::MLEstimator, x1::T, x2::T) where {T<:Integer}

    checkx1x2(x1, x2, design(estimator))
    return (x1 + x2)/samplesize(design(estimator), x1)

end
