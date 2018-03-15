[https://ci.appveyor.com/api/projects/status/github/kkunzmann/BinaryTwoStageDesigns]

# BinaryTwoStageDesigns

A [Julia](http://julialang.org) package implementing methods for planning and
evaluation of exact single-arm two-stage designs with binary
endpoint.
Documentation is available [here](https://imbi-heidelberg.github.io/BinaryTwoStageDesigns).

## Installation

Several key features depend on Integer Linear Programs and a suitable solver
for the [JuMP](https://github.com/JuliaOpt/JuMP.jl)-package must be installed
prior to using these.
It is strongly recommended to use Gurobi as solver 
[Gurobi](http://www.gurobi.com/index) and its Julia interface
[Gurobi.jl](https://github.com/JuliaOpt/Gurobi.jl).
Free academic licenses can be obtained
[here](http://www.gurobi.com/academia/for-universities).
Make sure Gurobi/Gurobi.jl/JuMP are working correctly before using any of the
functionality involving optimization.

Alternatively, the open source [Cbc solver](https://projects.coin-or.org/Cbc) can
be used for finding optimal designs ([JuMP integration](https://github.com/JuliaOpt/Cbc.jl)) 
and the [Ipopt solver](https://projects.coin-or.org/Ipopt) for minimizing the
quadratic objective of the optimal compatible estimator ([JuMP integration](https://github.com/JuliaOpt/Ipopt.jl))

The package is not yet registered and the current development verions can be
installed via:

    Julia> Pkg.clone("https://github.com/kkmann/BinaryTwoStageDesigns.git")

Note that the package is not yet final. 
A specific release can be installed manually by dowloading the release as .zip folder from 
https://github.com/imbi-heidelberg/BinaryTwoStageDesigns/releases and unpack to your
local .julia/v0.6 directory.



## Documentation

The documentation of the current development version is available
[here](https://kkmann.github.io/BinaryTwoStageDesigns/).
