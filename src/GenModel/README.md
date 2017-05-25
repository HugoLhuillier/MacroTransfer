## General model

To obtain the policy functions presented in the thesis,

```julia
p      = GenPFI.Param(size=20, ut = "CRRA", ɣ = 3., α = 0.9); # setup the model
pol    = GenPFI.ite(30,p;tol=1e-3)                            # iterate until convergence
```

### Model

Here, agents live for four periods: two of youth, and two of old ages. The stage game is the same as the one in the toy model. Due to the promising results obtained with policy function iteration in the toy model, and its easy extension to multi-dimensional problems, we use this method to solve numerically this model.

### Numerical challenges

**1)** *Larger state space and computational time*.

**2)** *Strategic interactions*. More periods means more strategic interactions. Mathematically, strategic interactions translate into the presence of the gradient of policy functions in the first order conditions. This, combined with occasionally binding constraints, can lead to discontinuous first order conditions.

As an example, Figure 1 plot the first order condition of the young 1, for a given state space. In the left panel, linear splines are used to interpolate the Markov strategies. In the right panel, cubic splines are used. In the second case, a non linear solver is able to find the root of the equation, while it cannot in the first one.

<p align="center">
  <b>Figure 1: First order condition, age 1</b>
  <br><br>
  <img src="https://github.com/HugoLhuillier/MacroTransfer/blob/master/output/figures/GenModel/NumProblems/FOC_smooth_nosmooth.png" alt="FOC 1" style="width: 400px;"/>
</p>

**3)** *Multiple-solutions* In the current version of the model, there is an infinity of possible allocation for the parent's wealth. For the oldest parent, age 4, we assume that all the wealth transmission occurs via intervivos transfer. For the age 3, we take the allocation in the middle of the solution set, as shown in Figure 2.

The old problems are multi-dimensional, non-linear root-finding problems, which are well known to be difficult to solve when little information is available on the functions. In particular, a solution may not exit, as shown in Figure 2.

<p align="center">
  <b>Figure 2: First order condition, age 3</b>
  <br><br>
  <img src="https://github.com/HugoLhuillier/MacroTransfer/blob/master/output/figures/GenModel/NumProblems/multiple_sol_3.png" alt="FOC 3" style="width: 400px;"/>
</p>

### Results

With the current calibration of the model, we obtain the following weatlh distributions.

<p align="center">
  <b>Figure 3: Wealth distributions</b>
  <br><br>
  <img src="https://github.com/HugoLhuillier/MacroTransfer/blob/master/output/figures/GenModel/wealth_histogram.png" alt="Wealth distributions" style="width: 300px;"/>
</p>
