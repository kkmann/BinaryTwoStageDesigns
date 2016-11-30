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
    conditionalProbabilityToReject, probabilityToReject

export optimalDesign, OptimalBinaryTwoStageDesign

export StageOneAdaptedOptimalBinaryTwoStageDesign,
    StageTwoAdaptedOptimalBinaryTwoStageDesign, adapt

include("BinaryTwoStageDesign.jl")
include("Parameters.jl")
include("util.jl")
include("OptimalBinaryTwoStageDesign.jl")
include("pointAlternative.jl")
include("adaptDesign.jl")

end # module
