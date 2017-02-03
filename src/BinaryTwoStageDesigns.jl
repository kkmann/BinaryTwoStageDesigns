__precompile__()

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


# include("Parameters.jl")
export Parameters, samplespace, maxsamplesize, efficacy, regularization,
    PointAlternative, null, alpha, alternative, beta,
    VagueAlternative, prior


# include("BinaryTwoStageDesign.jl")
export AbstractBinaryTwoStageDesign,
    BinaryTwoStageDesign,
        interimsamplesize, samplesize, criticalvalue, getRejectionBoundary,
        conditionalpower, power, test, pdf,
        simulate


# include("SampleSize.jl")
export SampleSize


# include("SampleSpace.jl")
export SampleSpace,
    SimpleSampleSpace, interimsamplesizerange, maxsamplesize, possible


abstract Regularization
type GroupSequential <: Regularization
end
type Unimodal <: Regularization
end
abstract Efficacy
type NoStoppingForEfficacy <: Efficacy
end
type StoppingForEfficacy <: Efficacy
end
export Regularization, GroupSequential, Unimodal,
    Efficacy, NoStoppingForEfficacy, StoppingForEfficacy

# include("SimpleMinimalExpectedSampleSize")
export SimpleMinimalExpectedSampleSize

# include optimal
export getoptimaldesign

export MinimalExpectedSampleSizePointAlternative,
    MinimalExpectedSampleSizeVagueAlternative


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

include("Parameters.jl")
include("BinaryTwoStageDesign.jl")
include("SampleSize.jl")
include("SampleSpace.jl")
include("Parameters/SimpleMinimalExpectedSampleSize.jl")

include("util.jl")
include("getoptimaldesign.jl")
# include("designs/pointAlternative.jl")
# include("adaptDesign.jl")
# include("estimate.jl")
# include("estimation/mle.jl")
# include("estimation/rb.jl")
# include("estimation/CompatibleEstimator.jl")
# include("ConfidenceInterval.jl")
# include("confidenceintervals/ClopperPearsonConfidenceInterval.jl")

end # module
