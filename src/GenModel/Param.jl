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
                  ρ           = 0.95,
                  lengthInc   = 2,
                  size        = 20,
                  ut          = "Log",
                  ɣ           = 1.0,
                  UI          = 2.0)

    # build the income stochastic process via Tauchen algorithm
    Y         = zeros(lengthInc, 4)
    П         = zeros(lengthInc, lengthInc,4)
    g(h::Int) = 1.5 * (0.3 * h^2 - 0.14*h^3/2 + 4.)
    s(h::Int) = h^2 - 5 * h + 7
    # for h in 1:4
    #   # mc      = QuantEcon.tauchen(lengthInc, ρ, h, g(h) - ρ/2 * g(h-1), 2)
    #   if h == 1
    #     mc      = QuantEcon.tauchen(lengthInc, ρ, s(h), g(h)+1, ceil(Int,h/2))
    #   else
    #     mc      = QuantEcon.tauchen(lengthInc, ρ, s(h), g(h) - ρ/2 * g(h-1), ceil(Int,h/2))
    #   end
    #   # Y[:,h]  = mc.state_values
    #   #if h == 4
    #     П[:,:,h]  = mc.p
    #   #end
    # end
    # Y[:,1]    = Y[:,2] / 2
    # Y[Y .< UI]= UI
    # Y    = Y
    # П         = [0.4 0.6; 0.2 0.8]
    П[:,:,1]    = [0.5 0.5; 0.4 0.6]
    П[:,:,2]    = [0.3 0.7; 0.2 0.8]
    П[:,:,3]    = [0.4 0.6; 0.1 0.9]
    П[:,:,4]    = [0.7 0.3; 0.3 0.7]
    for h in 1:4
      Y[1,h] = UI
      Y[2,h] = g(h)
    end

    a_grid    = linspace(1e-7, 4*maximum(Y), size)
    b_grid    = linspace(1e-7, 2*maximum(Y), size)

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
