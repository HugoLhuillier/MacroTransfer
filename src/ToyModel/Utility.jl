# returns the utility function
function pick_u(p::Param)
  gamma, ut      = p.gamma, p.ut

  if ut == "Log"
    du = (x) -> log(x)
  elseif ut == "CRRA"
    du = (x) -> (x^(1-gamma) - 1) / (1 - gamma)
  elseif ut == "Exp"
    du = (x) -> (1 - exp(-gamma * x)) / gamma
  elseif ut == "Quad"
    du = (x) -> x - gamma * x^2 / 2
  end

  return du
end

# returns the derivative of the utility function
function pick_du(p::Param)
  gamma, ut      = p.gamma, p.ut

  if ut == "Log"
    du = (x) -> 1 / x
  elseif ut == "CRRA"
    du = (x) -> x^(-gamma)
  elseif ut == "Exp"
    du = (x) -> Exp(-gamma * x)
  elseif ut == "Quad"
    du = (x) -> 1 - gamma * x
  end

  return du
end

# returns the inverse of the derivative.
function pick_idu(p::Param)
  if p.ut == "Log"
    idu = (x) -> 1 / x
  elseif p.ut == "CRRA"
    idu = (x) -> (1 / x)^(1 / p.gamma)
  elseif p.ut == "Exp"
    idu = (x) -> - log(x) / p.gamma
  elseif p.ut == "Quad"
    idu = (x) -> (1 - x) / p.gamma
  end
  return idu
end
