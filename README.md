# BinaryTwoStageDesigns

A [Julia](http://julialang.org) package implementing methods for planning and
evaluation of exact single-arm two-stage designs with binary
endpoint.
Documentation is available [here](https://imbi-heidelberg.github.io/BinaryTwoStageDesigns).

## Installation

Several key features depend on Integer Linear Programs and a suitable solver
for the [JuMP](https://github.com/JuliaOpt/JuMP.jl)-package must be installed
prior to using these.
The package is tested against the commercial solver
[Gurobi](http://www.gurobi.com/index) and its Julia interface
[Gurobi.jl](https://github.com/JuliaOpt/Gurobi.jl).
Free academic licenses for Gurobi can be obtained
[here](http://www.gurobi.com/academia/for-universities).
Make sure Gurobi/Gurobi.jl/JuMP are working correctly before using any of the
functionality involving optimization.

The package is not yet registered and the current development verions can be
installed via:

    Julia> Pkg.clone("https://github.com/imbi-heidelberg/BinaryTwoStageDesigns")

Note that these versions need not be stable (minor fixes and documentation changes). 
A specific release, which is just a tagged commit, must be checked out using git directly. 
I.e. switch to  ~/.julia/v0.6/BinaryTwoStageDesigns in your user folder and checkout the tag 
directly (git must be installed) by e.g.:

    git checkout v0.1.0

Note that this will cause a warning about a 'detached HEAD' state. This is
fine for usage of the specific version and only has implications for developers.



## Documentation

The documentation of the current development version is available
[here](https://imbi-heidelberg.github.io/BinaryTwoStageDesigns) and previous
versions can be found locally under /docs/build/index.html of the
respective release.
