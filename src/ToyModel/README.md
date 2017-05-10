# Toy model

To obtain the solution via the policy function algorithm, do

```julia
d = getToySol("PFI")

julia> Dict{String,Any} with 2 entries:
  "param" => # parameters
  "sol"   => # solutions
```

You can also modify the default parameters here, as well as the algorithm used -- either "VFI", "PFI" or "EGM".

### Model

In this version of the model, agents live only for two periods

  1. In the first period, newborns receive a given endowment, and a transfer from their parent.
  They then decide how much to consume and save out of it.
  2. In the second period, they consume and transfer a fraction of their wealth to their kid. They die at
  the end of the period.

This model is solved numerically both using value function iteration, policy function iteration, and the endogenous grid method of Carroll (2005). Additionally, when letting the model be deterministic, and some more assumptions (see pdf), we are able to derive closed form solutions for the Markov strategies, and the Markov perfect equilibrium. This allows us to compare the numerical results with the "true" solution.

### Results

Both PFI and EGM match the analytical solutions - see Figure 1. The number of iteration is fast, between 4 and 5 for both EGM and PFI depending on the initial guess.

<p align="center">
  <b>Figure 1: Numerical vs. analytical solutions</b>
  <br><br>
  <img src="https://github.com/HugoLhuillier/MacroTransfer/blob/master/output/figures/ToyModel/num_vs_analytical.png" alt="Numerical vs. analytical" style="width: 400px;"/>
</p>

On the contrary, VFI is much less accurate than the two other methods. This is due to

1. the presence of strategic interaction. Mathematically, this translates into the presence of the derivative of the policy function in the first order condition. When using EGM or PFI, we give this information to the algorithm, helping it to understand better the problem.  This is not the case when we are doing VFI.
1. the non-negativity constraint on transfers.

<p align="center">
  <b> Figure 2: Analytical vs. VFI </b>
  <br><br>
  <img src="https://github.com/HugoLhuillier/MacroTransfer/blob/master/output/figures/ToyModel/ana_vs_vfi.png" alt="Numerical vs. analytical" style="width: 400px;"/>
</p>

Finally, in agreement with the literature, we find that the fastest algorithm is EGM, followed by PFI. As expected, VFI is the least performant algorithm of those three.

<p align="center">
  <b>Figure 3: Performance</b>
  <br><br>
  <img src="https://github.com/HugoLhuillier/MacroTransfer/blob/master/output/figures/ToyModel/algo_perf.png" alt="Algo performance" style="width: 400px;"/>
</p>
