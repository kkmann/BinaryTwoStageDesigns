.. BinaryTwoStageDesigns documentation master file, created by
   sphinx-quickstart on Fri Jun 30 15:03:04 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Binary Two-Stage Designs
========================

A `Julia <http://julialang.org>`_ package implementing methods for planning and
evaluation of exact single-arm two-stage designs with binary endpoint.

.. module:: BinaryTwoStageDesigns
   :synopsis: Methods for planning and evaluation of exact single-arm two-stage
              designs with binary endpoint

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   quickstart.rst

`BinaryTwoStageDesigns <https://github.com/imbi-heidelberg/BinaryTwoStageDesigns>`_
provides functionality for working with exact single-arm two-stage designs for
clinical trials with binary endpoint.
Special focus lies on the computation of optimal two-stage designs extending
the ideas of [Simon1989]_ using the `Julia <http://julialang.org/>`_
programming language.
BinaryTwoStageDesigns relies heavily on the functionality of the
`JuMP <https://github.com/JuliaOpt/JuMP.jl>`_ package and is tested against the
commercial mixed integer programming solver
`Gurobi <http://www.gurobi.com/index>`_.
While JuMP supports other solvers the use of Gurobi is highly recommended.
Academic licenses for Gurobi can be obtained free of charge upon request.

.. 	[Simon1989] Simon, R. Optimal two-stage designs for phase II clinical trials.
  `Controlled Clinical Trials` 1989; 10, 1–10

Installation
------------

For those not familiar with Julia, `<http://julialang.org/learning/>`_ might be
a good starting point.
The language can be tried out without having to install any software via the
Julia containers provided by `<https://www.juliabox.org/>`_ for free
(Google account required for sign in).

Make sure Gurobi/`Gurobi.jl <https://github.com/JuliaOpt/Gurobi.jl>`/JuMP are
working correctly before using any of the functionality involving optimization.

The package is not yet registered and the current development verions can be
installed via

    Julia> Pgk.clone("https://github.com/imbi-heidelberg/BinaryTwoStageDesigns")

A specific release can then be installed by checking-out the respective
branch, e.g.

    Julia> Pkg.checkout("BinaryTwoStageDesigns", "0.1")

The documentation of the current development version is available under
`<https://imbi-heidelberg.github.io/BinaryTwoStageDesigns>` and previous versions
can be found under /docs/build/html/index.html of the respective release.

Quick-start
-----------

See :ref:`quickstart`.



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
