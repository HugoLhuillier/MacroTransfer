## General model

:rotating_light: :rotating_light: <b>Code in progress</b> :rotating_light: :rotating_light:

So far, one can use the model to compute manually the first iteration of the algorithm.

```julia
# initialize the model
p   = GenPFI.Param(size = 10)
U   = GenPFI.Utility(p)
pol = GenPFI.Policies(p)
# compute manually the first iteration, for the first three steps
GenPFI.old4!(p,U,pol)
GenPFI.young1!(p,U,pol)
GenPFI.old3!(p,U,pol)
GenPFI.young2!(p,U,pol)
```

### Model

Here, agents live for four periods: two of youth, and two of old ages. The stage game is the same as the one in the toy model. Due to the promising results obtained with policy function iteration in the toy model, and its easy extension to multi-dimensional problems, we use this method to solve numerically this model.

### Numerical challenges

**1)** *Large state and computational time*. For instance, for the age 2 problem, the size of the state space is NY x NY x NG x NG x NG, where NY denotes the number of grid points in the income grid, and NG the number of grid points in the asset and transfer grid. Table 1 reports the computational time for each problem, for different grid sizes.

<p align="center">
<b>TABLE 1: Computation time, in seconds </b>
<table class="tg">
  <tr>
    <th class="tg-031e"></th>
    <th class="tg-hgcj" colspan="3">Grid size</th>
  </tr>
  <tr>
    <td class="tg-031e"></td>
    <td class="tg-s6z2">5</td>
    <td class="tg-s6z2">10</td>
    <td class="tg-s6z2">20</td>
  </tr>
  <tr>
    <td class="tg-e3zv">Age 1</td>
    <td class="tg-s6z2">2.101</td>
    <td class="tg-s6z2">8.482</td>
    <td class="tg-s6z2">34.959</td>
  </tr>
  <tr>
    <td class="tg-e3zv">Age 3</td>
    <td class="tg-s6z2">0.089</td>
    <td class="tg-s6z2">0.443</td>
    <td class="tg-s6z2">1.237</td>
  </tr>
  <tr>
    <td class="tg-e3zv">Age 4</td>
    <td class="tg-s6z2">0.172</td>
    <td class="tg-s6z2">0.702</td>
    <td class="tg-s6z2">2.831</td>
  </tr>
</table>
</p>

**2)** *Strategic interactions*. More periods means more strategic interactions. Mathematically, strategic interactions translate into the presence of the gradient of policy functions in the first order conditions. This, combined with occasionally binding constraints, can lead to discontinuous first order conditions.

As an example, Figure 1 plot the first order condition of the young 1, for a given state space. In the left panel, linear splines are used to interpolate the Markov strategies. In the right panel, cubic splines are used. In the second case, a non linear solver is able to find the root of the equation, while it cannot in the first one.

<p align="center">
  <b>Figure 1: First order condition, age 1</b>
  <br><br>
  <img src="https://github.com/HugoLhuillier/MacroTransfer/blob/master/output/figures/GenModel/NumProblems/FOC_smooth_nosmooth.png" alt="FOC 1" style="width: 400px;"/>
</p>

**3)** *Multi-dimensionality* The old problems are multi-dimensional, non-linear root-finding problems, which are well known to be difficult to solve when little information is available on the functions. In particular, a solution may not exit, as shown in Figure 2.

<p align="center">
  <b>Figure 3: First order condition, age 3</b>
  <br><br>
  <img src="https://github.com/HugoLhuillier/MacroTransfer/blob/master/output/figures/GenModel/NumProblems/FOC_3_nosol.png" alt="FOC 3" style="width: 400px;"/>
</p>
