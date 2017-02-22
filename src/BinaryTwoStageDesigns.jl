__precompile__()

module BinaryTwoStageDesigns

using JuMP
using QuadGK

import Base.show
import Base.convert
import Distributions, Distributions.minimum, Distributions.maximum,
    Distributions.pdf, Distributions.cdf, Distributions.rand,
    Distributions.quantile, Distributions.mean, Distributions.var
import DataFrames, DataFrames.DataFrame
import DataArrays
import MathProgBase


# include("Parameters.jl")
export Parameters, samplespace, maxsamplesize, isgroupsequential, allowsstoppingforefficacy,
    NoParameters,
    PointAlternative, null, alpha, alternative, beta,
    VagueAlternative, prior


# include("BinaryTwoStageDesign.jl")
export AbstractBinaryTwoStageDesign,
    BinaryTwoStageDesign,
        parameters, interimsamplesize, samplesize, criticalvalue, power, test, pdf,
        simulate, support, ispossible, stoppingforfutility, score


# include("SampleSize.jl")
export SampleSize


# include("SampleSpace.jl")
export SampleSpace,
    SimpleSampleSpace, interimsamplesizerange, maxsamplesize, possible

# include("SimpleMinimalExpectedSampleSize.jl")
export SimpleMinimalExpectedSampleSize,
    IsGroupSequential, GroupSequential, NotGroupSequential,
    StoppingForEfficacy, AllowStoppingForEfficacy, NoStoppingForEfficacy,
    HasMonotoneConditionalPower, MonotoneConditionalPower, NoMonotoneConditionalPower,
    minconditionalpower, smoothness

export MinimalMinimalExpectedSampleSize

export LiuScore, ros, rup

export KunzmannScore, expectedcost, underpowerpenalty

# include ("optimal.jl")
export getoptimaldesign

# include ("adapt.jl")
export adapt

# include ("estimate.jl")
export BinaryTwoStageDesignEstimator,
    design, estimate, p, bias, rmse, incompatibleoutcomes,
    MaximumLikelihoodEstimator, # include("Estimators/mle.jl")
    RaoBlackwellizedEstimator, # inlcude("Estimators/rbe.jl")
    jeffreysprior, CompatibleEstimator

# include("")
export ConfidenceInterval,
    limits, confidence, design, coverage,
    ClopperPearsonConfidenceInterval, estimator

include("Parameters.jl")
include("BinaryTwoStageDesign.jl")
include("SampleSize.jl")
include("SampleSpace.jl")
include("Parameters/MinimalMinimalExpectedSampleSize.jl")
include("Parameters/SimpleMinimalExpectedSampleSize.jl")
include("Parameters/LiuScore.jl")
include("Parameters/KunzmannScore.jl")
include("util.jl")
include("getoptimaldesign.jl")
include("adapt.jl")
include("estimate.jl")
include("Estimators/mle.jl")
include("Estimators/rbe.jl")
include("Estimators/CompatibleEstimator.jl")
include("ConfidenceInterval.jl")
include("ConfidenceIntervals/ClopperPearsonConfidenceInterval.jl")

end # module
