# Master's thesis

[![Build Status](https://travis-ci.org/HugoLhuillier/MacroTransfer.svg?branch=master)](https://travis-ci.org/HugoLhuillier/MacroTransfer)
[![Build Status](https://img.shields.io/badge/Toy%20Model-Done-brightgreen.svg)](https://github.com/HugoLhuillier/MacroTransfer/tree/master/src/ToyModel)
[![Build Status](https://img.shields.io/badge/General%20Model-fine%20tuning-yellowgreen.svg)](https://github.com/HugoLhuillier/MacroTransfer/tree/master/src/GenModel)

This repository contains the code used for my Master's thesis at Sciences Po Paris.

In this repository, you'll find:

1. [src/ToyModel](https://github.com/HugoLhuillier/MacroTransfer/tree/master/src/ToyModel): uses value function iteration, policy function iteration, and the endogeneous grid method of Carroll (2005), to solves a simplified, two-periods version of the model. The solutions are also compared to the analytical ones, which can be derived in this case.
1. [src/GenModel](https://github.com/HugoLhuillier/MacroTransfer/tree/master/src/GenModel): uses policy function iteration to solve the larger model.
1. [output/env](https://github.com/HugoLhuillier/MacroTransfer/tree/master/output/env): contains the policy functions, the stationary distributions and the two-period panels used in the thesis.

:rotating_light:**Warning**:rotating_light:: this repository uses the `Interpolations.jl` package to interpolate multi-dimensional functions, which is itself based on `AxisAlgorithm.jl`. Its current version features a small bug, as mentioned [here](https://github.com/JuliaMath/Interpolations.jl/issues/151). The problem has been solved in [this pull request](https://github.com/timholy/AxisAlgorithms.jl/pull/7), but has not been merged yet. For the code in this repo to run, one thus needs to implement manually the solution given in the pull request.

## *Intergenerational transfers without commitment: a macroeconomic framework*
------------

We consider here an **heterogeneous agents model** with a lifecycle and overlapping generations. Each generation lives for several periods. In their first years of life, agents are considered as "youngs" - one can think of them as 20-40 year old individuals, who have finished their education, and have started to work. In the second stage of their life, they become "old". At this point, they all have a (unique) child, to whom they can transfer some wealth in every period. Transfers occur because olds are assumed to be altruistic towards their descendants.

In every periods, agents decide rationally how much to consume and save, while being prevented to borrow above some exogeneous limit, hence the emergence of heterogeneous agents. When old, they additionally choose the amount of wealth transferred to their child.

The key unit of interest in this model is the **family**. We depart from the unitary family model, and model it as a pair (young, old). Both members have their own preferences, and The decision-making process at the family level is modeled as a **non-cooperative, two-stage game**. In particular, parents cannot impose consumption and savings patterns to their child. Moreover, both members do not have access to **commitment devices**.

### A glimpse on results
#### Toy model

<p align="center">
  <b>Figure 1: Numerical vs. analytical solutions</b>
  <br><br>
  <img src="https://github.com/HugoLhuillier/MacroTransfer/blob/master/output/figures/ToyModel/num_vs_analytical.png" alt="Numerical vs. analytical" style="width: 400px;"/>
</p>

#### General model

<p align="center">
  <b>Figure 2: Wealth distributions</b>
  <br><br>
  <img src="https://github.com/HugoLhuillier/MacroTransfer/blob/master/output/figures/GenModel/wealth_histogram.png" alt="Wealth distributions" style="width: 300px;"/>
</p>
