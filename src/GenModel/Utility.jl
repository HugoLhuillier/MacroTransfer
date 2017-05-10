"""
```
Utility(p::Param)
```
Constructor of the `Utility` immutable, which contains the first, `du(x)`,
second, `du_p(x)`, and third, `du_s(x)`, order partial derivative of the policy functions.
"""
immutable Utility
  du::Function
  du_p::Function
  du_s::Function

  function Utility(p::Param)

    if p.ut == "Log"
      du   = function(x)
        if x > 0
          return 1 / x
        else
          return 1e9
        end
      end
      du_p = (x) -> - 1 / (x * x)
      du_s = (x) -> 2 / (x * x * x)
    elseif p.ut == "CRRA"
      du   = (x) -> x^(-p.ɣ)
      du_p = (x) -> -p.ɣ * x^(-1 - p.ɣ)
      du_s = (x) -> p.ɣ * (1 + p.ɣ) * x^(-2 - p.ɣ)
    elseif p.ut == "Exp"
      du   = (x) -> exp(-p.ɣ * x)
      du_p = (x) -> -p.ɣ * exp(-p.ɣ * x)
      du_s = (x) -> p.ɣ * p.ɣ * exp(-p.ɣ * x)
    else
      du   = (x) -> 1 - p.ɣ * x
      du_p = (x) -> - p.ɣ
      du_s = (x) -> 0
    end

    return new(du, du_p, du_s)
  end
end
