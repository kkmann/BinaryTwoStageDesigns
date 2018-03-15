using BinaryTwoStageDesigns
using Base.Test, QuadGK, Roots, Cbc, Ipopt
import Distributions

solver  = CbcSolver(seconds = 600)
qsolver = IpoptSolver()

include("test_Design.jl")

include("test_SampleSize.jl")

include("test_SampleSpace.jl")

include("test_MESS.jl")

include("test_MBESS.jl")

include("test_EB.jl")

include("test_Liu.jl")

include("test_optimization.jl")

include("test_estimators.jl")

include("test_ci.jl")
