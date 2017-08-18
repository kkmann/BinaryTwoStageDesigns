using BinaryTwoStageDesigns
using Base.Test
using QuadGK
using Gurobi

solver = GurobiSolver(
  MIPGap     = 10^(-4.0),
  TimeLimit  = 3600,
  OutputFlag = 0 
)

include("test_Design.jl")

include("test_SampleSize.jl")

include("test_SampleSpace.jl")

include("test_MESS.jl")

include("test_optimization.jl")

# include("test_estimators.jl")

# include("test_ci.jl")
