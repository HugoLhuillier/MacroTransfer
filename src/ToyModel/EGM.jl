#==========================#
# EGM functions
#==========================#

function young_E!(p::Param, _cy::Matrix, _co::Matrix)

  du, idu    = pick_du(p), pick_idu(p)
  C, bs      = zeros(size(_cy)), zeros(size(_cy))

  # Compute the right hand side of the FOC and invert it to obtain the consumption
  function foc_solver(a_idx::Integer, y_idx::Integer)
    T    = 0
    for k in 1:length(p.inc)
      T += p.Pi[y_idx, k] * du(_co[a_idx, k])
    end
    T   *= p.beta * p.R
    ce   = idu(T)
    b    = ce + p.a_grid[a_idx] - p.inc[y_idx]
    return ce, b
  end

  # Convert the function from the endo grid onto the initial grid
  function C_itp(b::Float64, y::Float64, yi::Integer)
    if b < bs[1, yi]
      cy  = b + y - p.a_grid[1]
    else
      cy  = linear_itp(b, bs[:,yi], C[:,yi])
    end
    return cy
  end

  for i in 1:length(p.a_grid)
    for j in 1:length(p.inc)
      C[i,j], bs[i,j] = foc_solver(i, j)
    end
  end

  for (k, b_k) in enumerate(p.b_grid)
    for (j, y_j) in enumerate(p.inc)
      _cy[k,j] = C_itp(b_k, y_j, j)
    end
  end

end

function old_E(p::Param, _cy::Matrix, _co::Matrix)

  du, idu     = pick_du(p), pick_idu(p)
  C, as, co   = zeros(size(_co)), zeros(size(_co)), zeros(size(_co))
  dcy         = der_pol(p, _cy)

  # compute the right hand side of the FOC and invert it
  function foc_solver(k::Integer, j::Integer)
    O   = p.alpha * dcy(p.b_grid[k], j) * du(_cy[k, j])
    ce  = idu(O)
    a   = (ce + p.b_grid[k]) / p.R
    return ce, a
  end

  # convert the function from the endo grid onto the initial grid
  function C_itp(a::Float64, yi::Integer)
    if a < as[1,yi]
      # corner solution
      c = p.R * a - p.b_grid[1]
    else
      # interior solution
      c = linear_itp(a, as[:,yi], C[:,yi])
    end
    return c
  end

  for k in 1:length(p.b_grid)
    for j in 1:length(p.inc)
      # solve for C and as
      C[k,j], as[k,j] = foc_solver(k, j)
    end
  end

  for (i, a_i) in enumerate(p.a_grid)
    for (j, y_j) in enumerate(p.inc)
      co[i,j]   = C_itp(a_i, j)
    end
  end

  return co

end

function init_E(p::Param)

  co      = zeros(length(p.a_grid), length(p.inc))
  cy      = zeros(length(p.b_grid), length(p.inc))
  for (y_idx, y) in enumerate(p.inc)
    co[:,y_idx] = p.a_grid / 2
    cy[:,y_idx] = p.beta * (p.b_grid + y) / (1 + p.beta)
  end

  return cy, co
end

function EGM(maxiter::Integer,
             p::Param;
             tol::Float64=1e-7,
             isCheck::Bool=false)

  cy, co_0  = init_E(p)

  for iter in 1:maxiter
    young_E!(p, cy, co_0)
    co_1    = old_E(p, cy, co_0)
    err     = maxabs(co_0 - co_1)

    if err < tol
      println("Found solutions after $iter iterations")
      # compute the two other Markov strategies
      a, b   = zeros(size(cy)), zeros(size(co_1))
      for (j, y) in enumerate(p.inc)
        a[:,j] = y + collect(p.b_grid) - cy[:,j]
        b[:,j] = p.R * collect(p.a_grid) - co_1[:,j]
      end

      return a, b, cy, co_1
      break
    elseif iter == maxiter
      error("EGM did not converge after $iter iterations")
    end

    if isCheck
      println("Iteration #$iter, error = ", err)
    end
    # update current guess
    co_0     = copy(co_1)
  end
end
