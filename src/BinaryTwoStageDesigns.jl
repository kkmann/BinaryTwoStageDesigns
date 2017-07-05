__precompile__()

module BinaryTwoStageDesigns

import JuMP
import QuadGK
import UnicodePlots
import Base.show, Base.print, Base.Multimedia.display
import Base.convert
import Distributions, Distributions.minimum, Distributions.maximum,
    Distributions.pdf, Distributions.cdf, Distributions.rand,
    Distributions.quantile, Distributions.mean, Distributions.var
import DataFrames, DataFrames.DataFrame
import DataArrays
import MathProgBase


# Parameters.jl
export Parameters, samplespace, maxsamplesize, isgroupsequential, allowsstoppingforefficacy,
    NoParameters,
    PointAlternative,
        null, alpha, alternative, beta


# BinaryTwoStageDesign.jl
export BinaryTwoStageDesign,
        parameters, interimsamplesize, samplesize, criticalvalue, power, test, pdf,
        simulate, support, ispossible, stoppingforfutility, score, expectedpower, jeffreysprior


# SampleSize.jl
export SampleSize


# SampleSpace.jl
export SampleSpace,
    SimpleSampleSpace, interimsamplesizerange, maxsamplesize, possible

export IPModel, extractsolution

# SimpleMinimalExpectedSampleSize.jl
export SimpleMinimalExpectedSampleSize,
    minconditionalpower, smoothness,
    getnvals, getcvals, isgroupsequential

# optimal.jl
export getoptimaldesign

# estimate.jl
export BinaryTwoStageDesignEstimator,
    design, estimate, p, bias, rmse, incompatibleoutcomes,
    MaximumLikelihoodEstimator, # include("Estimators/mle.jl")
    RaoBlackwellizedEstimator, # inlcude("Estimators/rbe.jl")
    CompatibleEstimator

#
export ConfidenceInterval,
    limits, confidence, design, coverage, estimator, findinconsistencies,
    meanwidth, meaninterval,
    ClopperPearsonConfidenceInterval,
    NaiveClopperPearsonConfidenceInterval,
    MinimumMeanWidthConfidenceInterval

include("util.jl")
include("Parameters.jl")
include("BinaryTwoStageDesign.jl")
include("SampleSize.jl")
include("SampleSpace.jl")
include("ipmodel.jl")
include("Parameters/SimpleMinimalExpectedSampleSize.jl")
include("getoptimaldesign.jl")
include("estimate.jl")
include("Estimators/mle.jl")
include("Estimators/rbe.jl")
include("Estimators/CompatibleEstimator.jl")
include("ConfidenceInterval.jl")
include("ConfidenceIntervals/ClopperPearsonConfidenceInterval.jl")
include("ConfidenceIntervals/NaiveClopperPearsonConfidenceInterval.jl")
include("ConfidenceIntervals/MinimumMeanWidthConfidenceInterval.jl")

end # module
