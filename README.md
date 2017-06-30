# BinaryTwoStageDesigns

A [Julia](http://julialang.org) package implementing methods for planning and
evaluation of exact single-arm two-stage designs with binary endpoint.

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
installed via

    Julia> Pgk.clone("https://github.com/imbi-heidelberg/BinaryTwoStageDesigns")

A specific release can then be installed by checking-out the respective
branch:

    Julia> Pkg.checkout("BinaryTwoStageDesigns", "0.1")

The documentation of the current development version is available under
https://imbi-heidelberg.github.io/BinaryTwoStageDesigns and previous versions
can be found under /docs/build/html/index.html of the respective release.
