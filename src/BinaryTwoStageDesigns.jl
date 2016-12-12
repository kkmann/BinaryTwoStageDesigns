module BinaryTwoStageDesigns

using JuMP

import Base.show
import Base.convert
import Distributions, Distributions.minimum, Distributions.maximum,
    Distributions.pdf, Distributions.cdf, Distributions.rand,
    Distributions.quantile, Distributions.mean, Distributions.var
import DataFrames, DataFrames.DataFrame
import DataArrays
import MathProgBase

export Parameters, PointAlternative,
    MinimalExpectedSampleSizePointAlternative, UncertainAlternative

export AbstractBinaryTwoStageDesign, BinaryTwoStageDesign,
    getInterimSampleSize, getSampleSize, SampleSize, getRejectionBoundary,
    conditionalProbabilityToReject, probabilityToReject, test, probability

export optimalDesign, OptimalBinaryTwoStageDesign

export StageOneAdaptedOptimalBinaryTwoStageDesign,
    StageTwoAdaptedOptimalBinaryTwoStageDesign, adapt

export BinaryTwoStageDesignEstimator, getDesign, estimate, p, bias, rmse

export MaximumLikelihoodEstimator

export RaoBlackwellizedEstimator

export CompatibleEstimator, jeffreysPrior

export BinaryTwoStageDesignConfidenceInterval, limits, getConfidence, getDesign,
    coverage

export ClopperPearsonConfidenceInterval

include("BinaryTwoStageDesign.jl")
include("Parameters.jl")
include("util.jl")
include("OptimalBinaryTwoStageDesign.jl")
include("designs/pointAlternative.jl")
include("adaptDesign.jl")
include("estimate.jl")
include("estimation/mle.jl")
include("estimation/rb.jl")
include("estimation/CompatibleEstimator.jl")
include("ConfidenceInterval.jl")
include("confidenceintervals/ClopperPearsonConfidenceInterval.jl")

end # module
