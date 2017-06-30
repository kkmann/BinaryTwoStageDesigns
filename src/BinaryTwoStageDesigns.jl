__precompile__()

module BinaryTwoStageDesigns

using JuMP
using QuadGK
import UnicodePlots
using Roots

import Base.show, Base.print
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
    PointAlternative, null, alpha, alternative, beta,
    VagueAlternative, prior, mcrv,
    simulate


# BinaryTwoStageDesign.jl
export AbstractBinaryTwoStageDesign,
    BinaryTwoStageDesign,
        parameters, interimsamplesize, samplesize, criticalvalue, power, test, pdf,
        simulate, support, ispossible, stoppingforfutility, score, expectedpower


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

export MinimalMinimalExpectedSampleSize

export LiuScore, ros, rup

export EB, expectedtransformedpower

export BESS

# optimal.jl
export getoptimaldesign

# estimate.jl
export BinaryTwoStageDesignEstimator,
    design, estimate, p, bias, rmse, incompatibleoutcomes,
    MaximumLikelihoodEstimator, # include("Estimators/mle.jl")
    RaoBlackwellizedEstimator, # inlcude("Estimators/rbe.jl")
    jeffreysprior, CompatibleEstimator

#
export ConfidenceInterval,
    limits, confidence, design, coverage, estimator, findinconsistencies,
    ClopperPearsonConfidenceInterval, NaiveClopperPearsonConfidenceInterval,
    MinimumMeanWidthConfidenceInterval

include("util.jl")
include("Parameters.jl")
include("BinaryTwoStageDesign.jl")
include("SampleSize.jl")
include("SampleSpace.jl")
include("ipmodel.jl")
include("Parameters/MinimalMinimalExpectedSampleSize.jl")
include("Parameters/SimpleMinimalExpectedSampleSize.jl")
include("Parameters/LiuScore.jl")
include("Parameters/BESS.jl")
include("Parameters/EB.jl")
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
