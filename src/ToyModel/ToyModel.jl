module ToyModel

using Interpolations, Roots, Optim
include("Param.jl")
include("Utility.jl")
include("Interp.jl")
include("PFI.jl")
include("EGM.jl")
include("VFI.jl")

export getSol

"""
```
getSol(algo::String;R=1.01, beta=0.98, alpha=0.8, Pi=[0. 1.;1. 0.], inc = [1., 6.], grid_size = 500, ut = "Log", gamma = 1.)
```

Solves the toy model for the parameters passed using

- value function iteration ("VFI")
- policy function iteration ("PFI")
- endogenous grid method ("EGM")

Returns dictionary containing

- param: parameters of the model
- sol: solutions of the model
"""
function getSol(algo::String;
                R         = 1.01,
                beta      = 0.98,
                alpha     = 0.8,
                Pi        = [0. 1.; 1. 0.],
                inc       = [1., 6.],
                grid_size = 100,
                ut        = "Log",
                gamma     = 1.0)


  p   = Param(R = R, beta = beta, alpha = alpha, Pi = Pi, inc = inc,
                grid_size = grid_size, ut = ut, gamma = gamma)
  if algo == "VFI"
    sol = VFI(100,p)
  elseif algo == "PFI"
    sol = PFI(100,p)
  elseif algo == "EGM"
    sol = EGM(100,p)
  else
    error("algo needs to be one of EGM, VFI or PFi")
  end
  return Dict("param" => p, "sol" => sol)
end

end
