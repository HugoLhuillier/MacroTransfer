#==========================#
# PFI functions
#==========================#

function young_P!(_B_c::Matrix,
                _A_c::Matrix,
                p::Param)

    R, beta, Pi    = p.R, p.beta, p.Pi
    inc, b_g       = p.inc, p.b_grid
    du             = pick_du(p)

    # Interpolate the transfer policy function & return FOC
    function foc_wrap(_y::Float64, _y_idx::Integer, _b::Float64)
      B_c = itp_pol(p, 'b', _B_c)
      foc = function(_x)
        e_du = 0
        for (yc_idx, yc) in enumerate(inc)
          if B_c[_x, yc_idx] > R * _x
            bb    = R * _x
          else
            bb    = B_c[_x, yc_idx]
          end
          if Pi[_y_idx, yc_idx] > 0
            e_du += du(R * _x - bb) * Pi[_y_idx, yc_idx]
          end
        end
        return du(_y + _b - _x) - beta * R * e_du
      end

      return foc
    end

    for (y_idx, y) in enumerate(inc)
      for (b_idx, b) in enumerate(b_g)

        # If currently poor, then kid will be rich. And vice versa
        foc   = foc_wrap(y, y_idx, b)
        if foc(0.) > 0
            # corner solution
            _A_c[b_idx, y_idx] = 0
        else
            # interior solution
            a_star             = fzero(foc, 0., y + b)
            _A_c[b_idx, y_idx] = a_star
        end
      end
    end
end

function old_P(_B_c::Matrix,
              _A_c::Matrix,
              p::Param)

  alpha, R       = p.alpha, p.R
  inc, a_g       = p.inc, p.a_grid
  du             = pick_du(p)
  B_c            = zeros(size(_B_c))

  # return the FOC as a function of the choice
  function foc_wrap(_a::Float64, _y::Float64, _yi::Integer)
    A_c       = itp_pol(p, 'a', _A_c)
    da = (_x) -> gradient(A_c, _x, _yi)
    foc       = function(_x::Float64)
      if A_c[_x, _yi] > _y + _x
        aa = _y + _x
      else
        aa = A_c[_x, _yi]
      end
      return du(R * _a - _x) - alpha * (1 - da(_x)[1]) * du(_y + _x - aa)
    end
    return foc
  end

  # solve for the optimal choice
  function foc_solver(_a::Float64, _y::Float64, _yi::Integer)
    foc   = foc_wrap(_a, _y, _yi)
    if foc(0.) > 0
      b_star = 0
    else
      b_star = fzero(foc, 0, R * _a)
    end
    return b_star
  end

  for (y_idx, y) in enumerate(inc)
    for (a_idx, a) in enumerate(a_g)
      B_c[a_idx, y_idx] = foc_solver(a, y, y_idx)
    end
  end

  return B_c
end

function init_P(p::Param)
  b_g, a_g, inc, R   = p.b_grid, p.a_grid, p.inc, p.R
  B_c       = zeros(length(a_g), length(inc))
  A_c       = zeros(length(b_g), length(inc))

  for (y_idx, y) in enumerate(inc)
    if y_idx == 1
      # Poor kid <-> rich parent
      B_c[:,y_idx] = a_g ./ 2
      A_c[:,y_idx] = (y .+ b_g) ./ (1 + R)
    else
      # Rich kid <-> poor parent
      B_c[:, y_idx] = 0.
      A_c[:, y_idx] = (y .+ b_g) ./ (1 + R)
    end
  end

  return A_c, B_c
end

function PFI(max_iter::Integer,
            p::Param;
            tol::Float64=1e-7,
            isCheck::Bool=false)

    A_c, B_c = init_P(p)

    for iter in 1:max_iter
      # update guess for young and old policy function, in that order
      young_P!(B_c, A_c, p)
      B_cn = old_P(B_c, A_c, p)
      err  = maxabs(B_cn - B_c)

      # check convergence
      if err < tol
        println("Found solutions after $iter iterations")

        cy, co = zeros(size(A_c)), zeros(size(B_cn))
        for (j, y) in enumerate(p.inc)
          cy[:,j] = y + collect(p.b_grid) - A_c[:,j]
          co[:,j] = p.R * collect(p.a_grid) - B_cn[:,j]
        end

        return A_c, B_cn, cy, co
        break
      elseif iter == max_iter
          error("No solution found after $iter iterations")
      end

      if isCheck
        println("Iteration # $iter. Error: ", err)
      end
      # update guess
      B_c = copy(B_cn)
    end
end
