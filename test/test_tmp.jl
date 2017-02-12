import BinaryTwoStageDesigns
import Gurobi
reload("BinaryTwoStageDesigns")

ss = BinaryTwoStageDesigns.SimpleSampleSpace(5:30, 75)
solver = Gurobi.GurobiSolver(
    IntFeasTol   = 1e-9,
    MIPGapAbs    = 1e-3,
    MIPGap       = 1e-4,
    Heuristics   = .25,
    NumericFocus = 3,
    Threads      = 2,
    MIPFocus     = 1,
    TimeLimit    = 900,
    OutputFlag   = 0
)
# params = BinaryTwoStageDesigns.MinimalMinimalExpectedSampleSize(
#     ss, .2, .4, .05, .2, .4
# )
# design1, res = BinaryTwoStageDesigns.getoptimaldesign(params, solver)
# println(convert(DataFrames.DataFrame, design1))
# params = BinaryTwoStageDesigns.SimpleMinimalExpectedSampleSize(
#     ss, .2, .4, .05, .2, .4, minconditionalpower = .7
# )
# design2, res = BinaryTwoStageDesigns.getoptimaldesign(params, solver)
# println(convert(DataFrames.DataFrame, design2))
# params = BinaryTwoStageDesigns.SimpleMinimalExpectedSampleSize(
#     ss, .2, .4, .05, .2, .4, minconditionalpower = .7, GROUPSEQUENTIAL = true
# )
# design3, res = BinaryTwoStageDesigns.getoptimaldesign(params, solver)
# println(convert(DataFrames.DataFrame, design3))
#
# params = BinaryTwoStageDesigns.SimpleMinimalExpectedSampleSize(
#     ss, .2, .4, .05, .2, .4, GROUPSEQUENTIAL = true, MONOTONECONDITIONALPOWER = false
# )
# design4, res = BinaryTwoStageDesigns.getoptimaldesign(params, solver)
# println(convert(DataFrames.DataFrame, design4))
ss = BinaryTwoStageDesigns.SimpleSampleSpace(1:10, 75)
params = BinaryTwoStageDesigns.LiuScore(
    ss, .2, .4, .05, .2, .5, .5
)
design3, res = BinaryTwoStageDesigns.getoptimaldesign(params, solver)
println(convert(DataFrames.DataFrame, design3))
