__precompile__()

module BinaryTwoStageDesigns

import JuMP
import QuadGK, Roots
import UnicodePlots
import Base.show, Base.print
import Base.convert, Base.writecsv
import Distributions, Distributions.minimum, Distributions.maximum,
    Distributions.pdf, Distributions.cdf, Distributions.rand,
    Distributions.quantile, Distributions.mean, Distributions.var
import DataFrames, DataFrames.DataFrame
import DataArrays
import MathProgBase
import CSV


# Parameters.jl
export 
  Parameters, 
    samplespace, maxsamplesize, isgroupsequential, allowsstoppingforefficacy, label, mcrv,
  NoParameters,
  PointAlternative,
    null, mtoer, alternative, mtter, mcrv,
  VagueAlternative,
    prior


# Design.jl
export 
  Design,
    parameters, interimsamplesize, samplesize, criticalvalue, power,
    expectedpower, test, pdf, simulate, support, ispossible,
    stoppingforfutility, stoppingforefficacy, score, expectedpower, jeffreysprior,
    save, writepropertiescsv


# SampleSize.jl
export 
  SampleSize


# SampleSpace.jl
export 
  SampleSpace,
    interimsamplesizerange, maxsamplesize, possible


# IPModel.jl
export 
  IPModel, 
    extractsolution

# MESS.jl
export 
  MESS,
    minconditionalpower, isgroupsequential

# LiuScore.jl
export 
  LiuScore,
    ros, rup

export 
  EB,
    expectedtransformedpower, expectedcost, expectedbenefit

export 
  MBESS, MBESS2

# optimal.jl
export optimaldesign

# estimate.jl
export 
  Estimator,
    design, estimate, pvalue, bias, rmse, mae, incompatibleoutcomes,
  MLEstimator, # include("Estimators/mle.jl")
  RBEstimator, # inlcude("Estimators/rbe.jl")
  OCEstimator

#
export 
  ConfidenceInterval,
    limits, confidence, design, coverage, estimator, findinconsistencies,
    meanwidth, meaninterval,
  ECPInterval,
  CPInterval,
  MMWInterval

export adaptstageone, adaptstagetwo

export adapt_stage_one, adapt_stage_two

include("util.jl")

include("Parameters.jl")

include("Design.jl")

include("SampleSize.jl")

include("SampleSpace.jl")

include("IPModel.jl")

include("Parameters/MESS.jl")
include("Parameters/LiuScore.jl")
include("Parameters/MBESS.jl")
include("Parameters/MBESS2.jl")
include("Parameters/EB.jl")

include("optimaldesign.jl")

include("Estimator.jl")
include("Estimators/MLEstimator.jl")
include("Estimators/RBEstimator.jl")
include("Estimators/OCEstimator.jl")

include("ConfidenceInterval.jl")
include("ConfidenceIntervals/ECPInterval.jl")
include("ConfidenceIntervals/CPInterval.jl")
include("ConfidenceIntervals/MMWInterval.jl")

include("adapt.jl")

end # module
