"""
```
Param(R=1.01, α=0.8, β=0.98, ρ=0.8, lengthInc=4, size=20, ut="Log", ɣ=1., UI=2.)
```

Constructor of the `Param` immutable. Output the `Param` immutable with fields

- `R::Float64`: Interest rate
- `α::Float64`: Altruism parameter
- `β:Float64`: Discount rate
- `П:Array{Float64}` : Transition matrix of the income Markov chain
- `Y::Array{Float64}` : Support of the income Markov chain
- `a_grid::LinSpace` : Asset grid
- `b_grid::LinSpace` : Transfer grid
- `ut::String` : Utility function. Needs to be one of Log, CRRA, Exp or Quad
- `ɣ::Float64` : Risk aversion parameter
"""
immutable Param
  R::Float64
  α::Float64
  β::Float64
  П::Array{Float64}
  Y::Array{Float64}
  a_grid::LinSpace
  b_grid::LinSpace
  ut::String
  ɣ::Float64

  function Param(;R           = 1.01,
                  α           = 0.8,
                  β           = 0.98,
                  ρ           = 0.8,
                  lengthInc   = 4,
                  size        = 20,
                  ut          = "Log",
                  ɣ           = 1.0,
                  UI          = 2.0)

    # build the income stochastic process via Tauchen algorithm
    Y         = zeros(lengthInc, 4)
    П         = zeros(lengthInc, lengthInc)
    g(h::Int) = (h^2 - 0.4 * h^3/2)
    for h in 2:4
      mc      = QuantEcon.tauchen(lengthInc, ρ, h, g(h) - ρ/2 * g(h-1), 2)
      Y[:,h]  = mc.state_values
      if h == 4
        П[:]  = mc.p
      end
    end
    Y[:,1]    = Y[:,2] / 2
    Y[Y .< UI]= UI
    Y    = Y

    a_grid    = linspace(1e-7, 2*maximum(Y), size)
    b_grid    = linspace(1e-7, maximum(Y), size)

    utility   = ["Log"; "CRRA"; "Exp"; "Quad"]
    if !(ut in utility)
        error("Utility needs to be one Log, CRRA, Exp or Quad")
    end
    if (ut == "Log") && (ɣ != 1.0)
        error("Log utility requires ɣ = 1.0")
    end

    return new(R, α, β, П, Y, a_grid, b_grid, ut, ɣ)
  end
end
