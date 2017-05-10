#==========================#
# VFI functions
#==========================#

function young_V!(_VO::Matrix, _A::Matrix, p::Param)
  VO = itp_pol(p, 'b', _VO)
  u  = pick_u(p)

  function obj_wrap(y::Float64, yi::Integer, b::Float64)
    f = function(x::Float64)
      VOp = 0
      for (_yi, _y) in enumerate(p.inc)
        VOp += p.Pi[yi, _yi] * VO[x, _yi]
      end
      return - u(y + b - x) - p.beta * VOp
    end
    return f
  end

  function max_solver(y::Float64, yi::Integer, b::Float64)
    f = obj_wrap(y, yi, b)
    return optimize(f, 0., y + b).minimizer
  end

  for (yi, y) in enumerate(p.inc)
    for (bi, b) in enumerate(p.b_grid)
      _A[bi, yi] = max_solver(y, yi, b)
    end
  end
end

function old_V!(_VO::Matrix, _A::Matrix, _B::Matrix, p::Param)
  A    = itp_pol(p, 'a', _A)
  VO_1 = zeros(size(_VO))
  u    = pick_u(p)

  function obj_wrap(y::Float64, yi::Integer, a::Float64)
    f = function(x::Float64)
      if y + x - A[x, yi] > 0
        return - u(p.R * a - x) - p.alpha * u(y + x - A[x, yi])
      else
        return Inf
      end
    end
    return f
  end

  function max_solver(y::Float64, yi::Integer, a::Float64)
    f = obj_wrap(y, yi, a)
    return optimize(f, 0., p.R * a).minimizer
  end

  for (yi, y) in enumerate(p.inc)
    for (ai, a) in enumerate(p.a_grid)
      _B[ai, yi]   = max_solver(y, yi, a)
      if y + _B[ai, yi] - A[_B[ai, yi], yi] < 0
        println("Here, error")
      end
      VO_1[ai, yi] = u(p.R * a - _B[ai, yi]) +
                    p.alpha * u(y + _B[ai, yi] - A[_B[ai, yi], yi])
    end
  end

  return VO_1
end

function init_V(p::Param)
  VO = zeros(length(p.a_grid), length(p.inc))
  A  = zeros(length(p.b_grid), length(p.inc))
  B  = zeros(length(p.a_grid), length(p.inc))

  for i in 1:length(p.inc)
    VO[:,i] = log.(p.R .* p.a_grid)
  end

  return VO, A, B
end

function VFI(maxIter::Integer, p::Param;
             tol::Float64 = 1e-6, isCheck::Bool = false)

  VO_0, A, B = init_V(p)
  for iter in 1:maxIter
    young_V!(VO_0, A, p)
    VO_1 = old_V!(VO_0, A, B, p)
    err  = maxabs(VO_1 - VO_0)

    if err < tol
      println("Found solutions after $iter iterations")
      return A, B
      break
    elseif iter == maxIter
      println("Did not find the solutions after $iter iterations.")
    end

    if isCheck
      println("Iteration #$iter: error = ", err)
    end

    VO_0 = copy(VO_1)
  end
end
