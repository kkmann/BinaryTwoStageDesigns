.. _quickstart:

Quick Start
===========

Make sure that Julia Gurobi/Gurobi.jl/JuMP and BinaryTwoStageDesigns are properly
installed on your system.

Setting up the Workspace
------------------------

.. code-block:: julia

  using BinaryTwoStageDesigns
  using Gurobi



Example problem
---------------

Assume that a pharmaceutical company wants to test a new anti-cancer agent in
a clinical phase II trial.
Tumor response after 30 days is adopted as surrogate endpoint for overall
survival.
A relatively large improvement under treatment over the  historical response
rate to treatment-as-usual (TAU) of :math:`p_0=0.2` is suspected.
Due to ethical concerns and randomization is considered infeasible and
the sponsor wants to reject the null hypothesis :math:`\mathcal{H}_0:p\leq p_0`
at a significance level of :math:`\alpha=5\%`.
At :math:`p_1=0.35` a power of at least :math:`80\%` is deemed sufficent.

Using these figures for planning a single-stage design (binomial test) would
result in a sample size of ????.
This has several disadvantages.
First of all, the discrete nature of the test statistic makes it impossible to
fully exhaust the type-I-error rate and is thus ineffective.
On the other hand, no safeguard agaist an actual lower response rate is possible
as the full number of subjects will always be included.

Two-stage designs can resolve these issues by allowing an interim analysis after
recruiting :math:`n_1` patients and deciding whether to continue to a second
stage or terminate the trial prematurely.
This idea was originally implemented by [Simon1989]_ using a simple
group-sequential structure (fixed size upon continuing).

Based on ideas presented in [Kunzmann2016]_ BinaryTwoStageDesigns
relaxes the constraint of a fixed sample size upon continuation and can be used
to derive designs with flexible stage-two sampel size which minimize the
expected sample size under a specific parameter value.

.. [Simon1989] Simon, R. Optimal two-stage designs for phase II clinical trials.
  `Controlled Clinical Trials` 1989; 10, 1–10

.. [Kunzmann2016] Kunzmann, K. and Kieser, K. Optimal adaptive two-stage designs for single-arm trials with binary endpoint.
  `https://arxiv.org/abs/1605.00249`


Finding Designs
---------------

For this example we will compute the design minimizing the expected sample size
under the alternative :math:`p_1=0.35` with the additional requirement that
at least 25 subjects should be enrolled if the trial does not stop early for
futility.
This constraint is sensible as it increases the reliability of the data on which
a potential subsequent phase III trial would be planned.

.. code-block:: julia

  SimpleSampleSpace()

todo


Exploring Design Properties
---------------------------

todo
