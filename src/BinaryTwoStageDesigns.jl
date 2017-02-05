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
export Parameters, samplespace, maxsamplesize, isgroupsequential, allowsstoppingforefficacy,
    NoParameters,
    PointAlternative, null, alpha, alternative, beta,
    VagueAlternative, prior


# include("BinaryTwoStageDesign.jl")
export AbstractBinaryTwoStageDesign,
    BinaryTwoStageDesign,
        interimsamplesize, samplesize, criticalvalue, getRejectionBoundary, power, test, pdf,
        simulate


# include("SampleSize.jl")
export SampleSize


# include("SampleSpace.jl")
export SampleSpace,
    SimpleSampleSpace, interimsamplesizerange, maxsamplesize, possible

# include("SimpleMinimalExpectedSampleSize")
export SimpleMinimalExpectedSampleSize,
    IsGroupSequential, GroupSequential, NotGroupSequential,
    StoppingForEfficacy, AllowStoppingForEfficacy, NoStoppingForEfficacy,
    minconditionalpower

# include optimal
export getoptimaldesign


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
