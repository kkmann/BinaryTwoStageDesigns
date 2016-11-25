module BinaryTwoStageDesigns

using JuMP

import Base.show
import Base.convert
import Distributions
import DataFrames
import MathProgBase

export PlanningAssumptions, optimalDesign

export BinaryTwoStageDesign,
    interimSampleSize, finalSampleSize, rejectionBoundary,
    conditionalProbabilityToReject, probabilityToReject,
    FinalSampleSize

export PointAlternative

include("BinaryTwoStageDesign.jl")
include("util.jl")
include("optimalDesign.jl")
include("pointAlternative.jl")

end # module
