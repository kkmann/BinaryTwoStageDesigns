module BinaryTwoStageDesigns

using JuMP

import Base.show
import Base.convert
import Distributions, Distributions.minimum, Distributions.maximum,
    Distributions.pdf, Distributions.cdf, Distributions.rand,
    Distributions.quantile, Distributions.mean, Distributions.var
import DataFrames
import MathProgBase

export PlanningAssumptions, optimalDesign

export BinaryTwoStageDesign
    getInterimSampleSize, getSampleSize, SampleSize, getRejectionBoundary,
    conditionalProbabilityToReject, probabilityToReject

export OptimalBinaryTwoStageDesign

export PointAlternative

include("BinaryTwoStageDesign.jl")
include("util.jl")
include("OptimalBinaryTwoStageDesign.jl")
include("pointAlternative.jl")

end # module
