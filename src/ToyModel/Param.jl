
"""
Type storing the information of the 2-period toy model.

##### Fields

- `R::Float64` : The interest rate plus 1 (strictly greater than 1)
- `beta::Float64` : Discount rate in (0, 1)
- `alpha:Float64` : Altruistic factor in (0, 1)
- `Pi::Matrix{Floa64}` : Transition matrix for newborn income
- `inc::Vector{Float64}` : Levels of income
- `b_grid::LinSpace{Float64}` : Transfer grid
- `a_grid::LinSpace{Float64}` : Asset grid
- `utility::String` : Utility function. Can choose between "Log", "CRRA", "Exp" and "Quad"
- `gamma::Float64` : Coefficient of risk aversion.
"""
type Param
  R::Float64
  beta::Float64
  alpha::Float64
  Pi::Array{Float64}
  inc::Array{Float64}
  b_grid::AbstractArray
  a_grid::AbstractArray
  ut::String
  gamma::Float64
end

function Param(;R         = 1.01,
                beta      = 0.98,
                alpha     = 0.8,
                Pi        =  [0. 1.;
                              1. 0.],
                inc       = [1., 6.],
                max_a     = Void,
                max_b     = Void,
                grid_size = 100,
                ut        = "Log",
                gamma     = 1.0)

  utility = ["Log"; "CRRA"; "Exp"; "Quad"]
  if !(ut in utility)
      error("Utility needs to be one of Log, CRRA, Exp or Quad")
  end
  if (ut == "Log") && (gamma != 1.0)
      error("Log utility requires gamma = 1")
  end

  if max_a == Void
    max_a  = 2 * maximum(inc)
  end
  if max_b == Void
    max_b  = maximum(inc)
  end
  b_grid = linspace(1e-7, max_b, grid_size)
  a_grid = linspace(1e-7, max_a, grid_size)

  return Param(R, beta, alpha, Pi, inc, b_grid, a_grid, ut, gamma)
end
