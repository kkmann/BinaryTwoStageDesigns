.. BinaryTwoStageDesigns documentation master file, created by
   sphinx-quickstart on Sun Dec 11 12:07:54 2016.

Binary Two-Stage Designs
========================

.. module:: AdaptiveDesigns
   :synopsis: Methods for Adaptive Clinical Trial Designs in Julia

`BinaryTwoStageDesigns <https://github.com/kkmann/BinaryTwoStageDesigns>`_
provides functionality for working with single arm two-stage designs for
clinical trials with binary endpoint.
Special focus lies on the computation of optimal two-stage designs extending
the ideas of [Simon1989]_ using the `Julia <http://julialang.org/>`_
programming language.
The **AdaptiveDesigns** package relies heavily on the functionality of the
`JuMP <https://github.com/JuliaOpt/JuMP.jl>`_ package and the mixed integer
programming solver `Gurobi <http://www.gurobi.com/index>`_.
While JuMP supports other solvers the use of Gurobi is highly recommended.
Academic licenses for Gurobi can be obtained free of charge upon request.

.. 	[Simon1989] Simon, R. Optimal two-stage designs for phase II clinical trials.
	`Controlled Clinical Trials` 1989; 10, 1–10

Installation
------------

A working Julia and Gurobi installation are required.
Simply add the package repository by::

	Pkg.clone("https://github.com/kkmann/BinaryTwoStageDesigns")

For those not familiar with Julia, `<http://julialang.org/learning/>`_ might be
a good starting point.
The language can be tried out without having to install any software via the
Julia containers provided by `<https://www.juliabox.org/>`_ for free
(Google account required for sign in).
When you have a working Julia setup running and you want to dive right into the
package usage see :ref:`quickstart`.



.. toctree::
   :maxdepth: 2
   :caption: Contents:

   quickstart.rst


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
