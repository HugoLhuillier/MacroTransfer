# Master's thesis

[![Build Status](https://travis-ci.org/HugoLhuillier/MacroTransfer.svg?branch=master)](https://travis-ci.org/HugoLhuillier/MacroTransfer)
[![Build Status](https://img.shields.io/badge/Toy%20Model-Done-brightgreen.svg)](https://github.com/HugoLhuillier/MacroTransfer/tree/master/src/ToyModel)
[![Build Status](https://img.shields.io/badge/General%20Model-In%20progress-orange.svg)](https://github.com/HugoLhuillier/MacroTransfer/tree/master/src/GenModel)

This repository contains the code used for my Master's thesis at Sciences Po Paris.

This repository is structured in the following way:

1. [src/ToyModel](https://github.com/HugoLhuillier/MacroTransfer/tree/master/src/ToyModel): uses value function iteration, policy function iteration, and the endogeneous grid method of Carroll (2005), to solves a simplified, two-periods version of the model. The solutions are also compared to the analytical ones, which can be derived in this case.
1. [src/GenModel](https://github.com/HugoLhuillier/MacroTransfer/tree/master/src/GenModel): uses policy function iteration to solve the larger model. This part is still work in progress.

:rotating_light:**Warning**:rotating_light:: this repository uses the `Interpolations.jl` package to interpolate multi-dimensional functions, which is itself based on `AxisAlgorithm.jl`. Its current version features a small bug, as mentioned [here](https://github.com/JuliaMath/Interpolations.jl/issues/151). The problem has been solved in [this pull request](https://github.com/timholy/AxisAlgorithms.jl/pull/7), but has not been merged yet. For the code in this repo to run, one thus needs to implement manually the solution given in the pull request.

## *Intergenerational transfers without commitment: a macroeconomic framework*
------------

We consider here an **heterogeneous agents model** with overlapping generations. Each generation lives for several periods. In their first years of life, agents are considered as "youngs" - one can think of them as 20-40 year old individuals, who have finished their education, and have started to work. In the second stage of their life, they become "old". At this point, they all have a (unique) child, to whom they can transfer some wealth in every period. Transfers occur because olds are assumed to be altruistic towards their descendants.

In every periods, agents decide rationally how much to consume and save, while being prevented to borrow above some exogeneous limit, hence the emergence of heterogeneous agents. When old, they additionally choose the amount of wealth transferred to their child.

A key assumption of this model is the **impossibility for parents to commit** to future transfers. This assumption is in line with the imperfect risk sharing observed empirically within families. As such, parents cannot directly control the consumption and savings of their descendants. Instead, they can only affect their behaviors via their own decisions. This leads to strategic interactions. Moreover, a direct consequence of this assumption is that transfers cannot be negative.

As it is well known in the literature, this kind of models are two complex to be solved analytically - except in some special cases, and we rely on numerical methods to find the optimal behaviors of agents. Strategic interactions complicate further the task, and their implications are discussed in length in the paper.

<p align="center">
<b>TABLE 1: Data vs. heterogenous agent models, De Nardi (2016)</b>
<table class="tg">
  <tr>
    <th class="tg-031e"></th>
    <th class="tg-s6z2"></th>
    <th class="tg-s6z2" colspan="5">Percentage wealth in the top</th>
  </tr>
  <tr>
    <td class="tg-031e"></td>
    <td class="tg-s6z2">Wealth<br>gini</td>
    <td class="tg-s6z2">1%</td>
    <td class="tg-s6z2">5%</td>
    <td class="tg-baqh">20%</td>
    <td class="tg-baqh">40%</td>
    <td class="tg-baqh">60%</td>
  </tr>
  <tr>
    <td class="tg-031e">U.S. data, 1989</td>
    <td class="tg-s6z2">0.78</td>
    <td class="tg-s6z2">29%</td>
    <td class="tg-s6z2">53%</td>
    <td class="tg-baqh">80%</td>
    <td class="tg-baqh">93%</td>
    <td class="tg-baqh">98%</td>
  </tr>
  <tr>
    <td class="tg-031e">No intergenerational links</td>
    <td class="tg-s6z2">0.67</td>
    <td class="tg-s6z2">7%</td>
    <td class="tg-s6z2">27%</td>
    <td class="tg-s6z2">69%</td>
    <td class="tg-s6z2">90%</td>
    <td class="tg-s6z2">98%</td>
  </tr>
  <tr>
    <td class="tg-031e">Parent's bequest motive</td>
    <td class="tg-s6z2">0.74</td>
    <td class="tg-s6z2">14%</td>
    <td class="tg-s6z2">37%</td>
    <td class="tg-s6z2">76%</td>
    <td class="tg-s6z2">90%</td>
    <td class="tg-s6z2">98%</td>
  </tr>
  <tr>
    <td class="tg-031e">Parent's bequest motive +<br>productivity inheritance</td>
    <td class="tg-s6z2">0.76</td>
    <td class="tg-s6z2">18%</td>
    <td class="tg-s6z2">42%</td>
    <td class="tg-s6z2">79%</td>
    <td class="tg-s6z2">95%</td>
    <td class="tg-s6z2">100%</td>
  </tr>
</table>
</p>


**The motivation** to this model is twofold

1. Since their apparition in economics, heterogeneous agents models have struggled to match the wealth distribution, as shown in Table 1. Including intergenerational transfers, in addition to bequests, can be one of the way to improve them.
1. Strategic interactions are very rarely included in macroeconomic models, despite their empirical relevance. The family is one of the most salient place in which these kind of behaviors occur.

### A glimpse on results
#### Toy model

<p align="center">
  <b>Figure 1: Numerical vs. analytical solutions</b>
  <br><br>
  <img src="https://github.com/HugoLhuillier/MacroTransfer/blob/master/output/figures/ToyModel/num_vs_analytical.png" alt="Numerical vs. analytical" style="width: 400px;"/>
</p>
