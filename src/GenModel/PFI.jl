module GenPFI

using Interpolations, ForwardDiff, NLsolve
using QuantEcon
using Logging
Logging.configure(level=INFO)
include("Param.jl")
include("Policies.jl")
include("Interp.jl")
include("Utility.jl")

#==========================#
# PFI functions
#==========================#

function young1!(p::Param, U::Utility, pol::Policies)
  # interpolate the future and other agents markov strategies
  B4, A4 = interp(p, pol.B4, "B4"), interp(p, pol.A4, "A4")
  B3, A3 = interp(p, pol.B3, "B3"), interp(p, pol.A3, "A3")
  A2     = interp(p, pol.A2, "A2")

  # FOC, to be solved for
  function f!(x, fvec, a_fp::Float64, b::Float64, y_c::Float64, yi_c::Int, yi_p::Int)
    # ensures that one cannot consume more than what one's have
    # less costly than putting an upper bound on the solver
    if x[1] > y_c + b
      fvec[1] = 1e9
    else
      Eu2 = 0
      Eu3 = 0

      # compute the t = 2 expectation
      for yi_fc in 1:length(p.Y[:,2])
        for yi_fp in 1:length(p.Y[:,4])
          _a4  = A4[yi_fp,yi_fc][a_fp, x[1]]
          _b4  = B4[yi_fp, yi_fc][a_fp, x[1]]
          _a2  = A2[yi_fc][x[1], _a4, _b4]

          # NOTE: gradient() do not access type ForwardDiff.Dual
          if typeof(x[1]) <: ForwardDiff.Dual{1,Float64}
            _da  = gradient(A4[yi_fp,yi_fc], a_fp, x[1].value)[2]
            _db  = gradient(B4[yi_fp,yi_fc], a_fp, x[1].value)[2]
          else
            _da  = gradient(A4[yi_fp,yi_fc], a_fp, x[1])[2]
            _db  = gradient(B4[yi_fp,yi_fc], a_fp, x[1])[2]
          end

          Eu2 += p.П[yi_c, yi_fc] * p.П[yi_p, yi_fp] *
                U.du(p.R * x[1] + p.Y[yi_fc, 2] + _b4 - _a2) * (p.R + _db)

          # compute the t = 3 expectation
          for yi_ffp in 1:length(p.Y[:,3])
            for yi_ffc in 1:length(p.Y[:,1])
              _a3 = A3[yi_ffp,yi_ffc][_a2 + _a4]
              _b3 = B3[yi_ffp,yi_ffc][_a2 + _a4]
              Eu3 += p.П[yi_c, yi_fc] * p.П[yi_fc, yi_ffp] * p.П[yi_fc, yi_ffc] * p.П[yi_p, yi_fp] *
                    U.du(p.R * (_a2 + _a4) + p.Y[yi_ffp,3] - _a3 - _b3)
            end
          end
          Eu3 *= _da
        end
      end

      fvec[1] = U.du(y_c + b - x[1]) - p.β * (Eu2 + p.R * p.β * Eu3)
    end
  end

  function f(x, a_fp::Float64, b::Float64, y_c::Float64, yi_c::Int, yi_p::Int)
    Eu2 = 0
    Eu3 = 0

    # compute the t = 2 expectation
    for yi_fc in 1:length(p.Y[:,2])
      for yi_fp in 1:length(p.Y[:,4])
        _a4  = A4[yi_fp,yi_fc][a_fp, x]
        _b4  = B4[yi_fp, yi_fc][a_fp, x]
        _a2  = A2[yi_fc][x, _a4, _b4]

        _da  = gradient(A4[yi_fp,yi_fc], a_fp, x)[2]
        _db  = gradient(B4[yi_fp,yi_fc], a_fp, x)[2]

        Eu2 += p.П[yi_c, yi_fc] * p.П[yi_p, yi_fp] *
              U.du(p.R * x + p.Y[yi_fc, 2] + _b4 - _a2) * (p.R + _db)

        # compute the t = 3 expectation
        for yi_ffp in 1:length(p.Y[:,3])
          for yi_ffc in 1:length(p.Y[:,1])
            _a3 = A3[yi_ffp,yi_ffc][_a2 + _a4]
            _b3 = B3[yi_ffp,yi_ffc][_a2 + _a4]
            Eu3 += p.П[yi_c, yi_fc] * p.П[yi_fc, yi_ffp] * p.П[yi_fc, yi_ffc] * p.П[yi_p, yi_fp] *
                  U.du(p.R * (_a2 + _a4) + p.Y[yi_ffp,3] - _a3 - _b3)
          end
        end
        Eu3 *= _da
      end
    end

    return U.du(y_c + b - x) - p.β * (Eu2 + p.R * p.β * Eu3)
  end

  function solver(a_fp::Float64, b::Float64, y_c::Float64, yi_c::Int, yi_p::Int)
    # find the roots of the FOC
    # with non-monotonic function, mcpsolve will find that the boundary condition binds
    # only if start from 0. however, it also finds sometimes that it binds, even though
    # it does not (e.g. if f(0) > 0, f'(0) > 0 but f(.) = 0 exists due to non-monotonicity)
    # to circumvent this: (i) check that there is f(.) > 0 everywhere
    # (ii.a) if so, then corner binds
    # (ii.b) else, looking for an interior solution, with multiple starting values if necessary
    # NOTE: fzero from Roots would allow us to not rely on multiple starting values. it is however
    # considerably slower than nlsolve
    debug("New solver, ($a_fp, $b, $y_c, $yi_c, $yi_p)")

    X = linspace(0., b+y_c-1e-1,100)
    F = ones(length(X))
    F = f.(X, a_fp,b,y_c,yi_c,yi_p)
    # if any(x -> x == true, F .< 0)
    #   debug("Interior")
    #   r   = fzero(x -> f(x,a_fp,b,y_c,yi_c,yi_p),b/2,[0.;b + y_c])
    #   if length(r) > 1
    #     debug("Multiple solution")
    #     _a1 = r[end]
    #   else
    #     _a1 = r[1]
    #   end
    # else
    #   debug("Corner")
    #   _a1 = 0.
    # end
    # return _a1

    if any(x -> x == true, F .< 0)
      debug("Interior")
      # r = NLsolve.nlsolve((x, fvec) -> f!(x, fvec, a_fp, b, y_c, yi_c, yi_p), [0.], [Inf],
      #            [(b + y_c)/(1+p.R)], reformulation = :smooth, iterations = 150, autodiff = true)
      r = NLsolve.nlsolve((x, fvec) -> f!(x, fvec, a_fp, b, y_c, yi_c, yi_p),
                  [(b + y_c)/(1+p.R)], iterations = 150)
      if !converged(r)
        debug("Not converged by middle")
        r = NLsolve.nlsolve((x, fvec) -> f!(x, fvec, a_fp, b, y_c, yi_c, yi_p),
                      [0.], iterations = 100)
        if !converged(r)
          debug("Not converged by 0")
          r = NLsolve.nlsolve((x, fvec) -> f!(x, fvec, a_fp, b, y_c, yi_c, yi_p),
                        [b + y_c - 1e-2], iterations = 100)
          if !converged(r)
            warn("Did not converged, ($a_fp, $b, $y_c, $yi_c, $yi_p)")
          end
        end
      end
      _a1 = r.zero[1]
    else
      debug("Corner")
      _a1 = 0.
    end

    return _a1
  end

  for (ai_fp, a_fp) in enumerate(p.a_grid)
    for (bi, b) in enumerate(p.b_grid)
      for (yi_c, y_c) in enumerate(p.Y[:,1])
        for (yi_p, y_p) in enumerate(p.Y[:,3])
          _A1 = solver(a_fp, b, y_c, yi_c, yi_p)
          # _A1 may be negative due to the cubic interpolation
          pol.A1[ai_fp, bi, yi_c, yi_p] = _A1 >= 0 ? _A1 : 0.
        end
      end
    end
  end
  info("Young1! ended")

end

function young2!(p::Param, U::Utility, pol::Policies)
  A3, B3  = interp(p, pol.A3, "A3"), interp(p, pol.B3, "B3")

  function f!(x, fvec, a_c::Float64, a_fp::Float64, b::Float64, y_c::Float64, yi_c::Int)
    # ensures that one cannot consume more than what one has
    if x[1] > p.R * a_c + b + y_c
      fvec[1] = 1e9
    else
      Eu3 = 0
      # compute the t = 2 expectation
      for yi_fp in 1:length(p.Y[:,3])
        for yi_fc in 1:length(p.Y[:,1])
          _b3 = B3[yi_fp,yi_fc][x[1] + a_fp]
          _a3 = A3[yi_fp,yi_fc][x[1] + a_fp]
          Eu3 += U.du(p.R * (x[1] + a_fp) + p.Y[yi_fp,3] - _b3 - _a3) * p.П[yi_c,yi_fp] * p.П[yi_c,yi_fc]
        end
      end
      fvec[1] = U.du(p.R * a_c + b + y_c - x[1]) - p.R * p.β * Eu3
    end
  end

  function f(x, a_c::Float64, a_fp::Float64, b::Float64, y_c::Float64, yi_c::Int)
    # ensures that one cannot consume more than what one has
    if x > p.R * a_c + b + y_c
      return 1e9
    else
      Eu3 = 0
      # compute the t = 2 expectation
      for yi_fp in 1:length(p.Y[:,3])
        for yi_fc in 1:length(p.Y[:,1])
          _b3 = B3[yi_fp,yi_fc][x + a_fp]
          _a3 = A3[yi_fp,yi_fc][x + a_fp]
          Eu3 += U.du(p.R * (x + a_fp) + p.Y[yi_fp,3] - _b3 - _a3) * p.П[yi_c,yi_fp] * p.П[yi_c,yi_fc]
        end
      end
      return U.du(p.R * a_c + b + y_c - x) - p.R * p.β * Eu3
    end
  end

  function solver(a_c::Float64, a_fp::Float64, b::Float64, y_c::Float64, yi_c::Int)
    # find the root of the FOC. use same method as young!1, see above
    X = linspace(0., p.R * a_c + b + y_c - 1e-1,100)
    F = ones(length(X))
    F = f.(X, a_c,a_fp,b,y_c,yi_c)
    if any(x -> x == true, F .< 0)
      debug("Interior")
      # r = NLsolve.mcpsolve((x, fvec) -> f!(x, fvec, a_c, a_fp, b, y_c, yi_c), [0.], [Inf],
      #              [(p.R * a_c + b + y_c) / (1 + p.R)], reformulation = :smooth,
      #              iterations = 150, autodiff = true)
      r = NLsolve.nlsolve((x, fvec) -> f!(x, fvec, a_c, a_fp, b, y_c, yi_c),
                   [(p.R * a_c + b + y_c) / (1 + p.R)], iterations = 150)
      if !converged(r)
        debug("Not converged by middle")
        r = NLsolve.nlsolve((x, fvec) -> f!(x, fvec, a_c, a_fp, b, y_c, yi_c),
                     [0.], iterations = 150)
        if !converged(r)
          debug("Not converged by 0")
          # r = NLsolve.mcpsolve((x, fvec) -> f!(x, fvec, a_c, a_fp, b, y_c, yi_c), [0.], [Inf],
          #               [p.R * a_c + b + y_c - 1e-2], reformulation = :smooth, iterations = 100, autodiff = true)
          r = NLsolve.nlsolve((x, fvec) -> f!(x, fvec, a_c, a_fp, b, y_c, yi_c),
                        [p.R * a_c + b + y_c - 1e-2], iterations = 150)
          if !converged(r)
            warn("Did not converged, ($a_fp, $b, $y_c, $yi_c, $yi_p)")
          end
        end
      end
      _a1 = r.zero[1]
    else
      debug("Corner")
      _a1 = 0.
    end
    return _a1
  end

  err = 0
  for (ai_c, a_c) in enumerate(p.a_grid)
    for (ai_fp, a_fp) in enumerate(p.a_grid)
      for (bi, b) in enumerate(p.b_grid)
        for (yi_c, y_c) in enumerate(p.Y[:,2])
            _A2 = solver(a_c,a_fp,b,y_c,yi_c)
            _A2 = _A2 >= 0. ? _A2 : 0.
            # update the distance between the two policy functions
            err = abs(pol.A2[ai_c,ai_fp,bi,yi_c] - _A2) > err ?
                  abs(pol.A2[ai_c,ai_fp,bi,yi_c] - _A2) : err
            pol.A2[ai_c,ai_fp,bi,yi_c] = _A2
        end
      end
    end
  end
  info(" Young2! ended")

  return err
end

function old3!(p::Param, U::Utility, pol::Policies)
  # markov strategies to be interpolated
  A4, B4  = interp(p, pol.A4, "A4"), interp(p, pol.B4, "B4")
  A1      = interp(p, pol.A1, "A1")
  A2      = interp(p, pol.A2, "A2")

  function f!(x, fvec, a_p::Float64, y_p::Float64, y_c::Float64, yi_p::Int, yi_c::Int)
    if x[1] + x[2] > p.R * a_p + y_p
      fvec[1] = 1e9
      fvec[2] = 1e9
    else
      _a1          = A1[yi_c, yi_p][x[1], x[2]]
      if typeof(x[1]) <: ForwardDiff.Dual{2,Float64}
        _da_a, _da_b = gradient(A1[yi_c, yi_p], x[1].value, x[2].value)
      else
        _da_a, _da_b = gradient(A1[yi_c, yi_p], x[1], x[2])
      end
      # compute the t = 4 expectations
      Eu4, Eu2 = 0., 0.
      for yi_fp in 1:length(p.Y[:,4])
        for yi_fc in 1:length(p.Y[:,2])
          _a4  = A4[yi_fp,yi_fc][x[1],_a1]
          _b4  = B4[yi_fp,yi_fc][x[1],_a1]
          _a2  = A2[yi_fc][_a1,_a4,_b4]

          Eu4 += U.du(p.R * x[1] + p.Y[yi_fp,4] - _b4 - _a4) * p.П[yi_p,yi_fp] * p.П[yi_c,yi_fc]
          Eu2 += U.du(p.R * _a1 + p.Y[yi_fc,2] + _b4 - _a2) * p.П[yi_p,yi_fp] * p.П[yi_c,yi_fc]
        end
      end
      fvec[1] = U.du(p.R * a_p + y_p - x[2] - x[1]) + p.α * _da_a * U.du(y_c + x[2] - _a1) -
              p.R * p.β * (Eu4 + p.α * _da_a * Eu2)
      fvec[2] = U.du(p.R * a_p + y_p - x[2] - x[1]) - p.α * (1 - _da_b) * U.du(y_c + x[2] - _a1) -
              p.α * p.R * p.β * _da_b * Eu2
    end
  end

  function solver(a_p::Float64, y_p::Float64, y_c::Float64, yi_p::Int, yi_c::Int, ai_p::Int)
    # find the roots of the system of equations. try to guess in a smart way the
    # initial value
    c0 = (p.R * a_p + y_p) / 7
    a0 = min(p.R * a_p + y_p - c0, y_p + a_p / (1 + p.R))
    b0 = max(0., p.R * a_p + y_p - c0 - a0)
    r = NLsolve.mcpsolve((x, fvec) -> f!(x, fvec, a_p, y_p, y_c, yi_p, yi_c), [0., 0.], [Inf, Inf],
                 [a0, b0], reformulation = :smooth, iterations = 150, autodiff = true)
    # in most of the cases, the solution does exist, but mcpsolve did not find it => multi-start approach
    i = 1
    while (!converged(r)) && (i <= 10)
      if i == 1
        # NOTE: to find the roots, we use instead nlsolve, as it works better than mcpsolve
        # also, note that we can be sure that the solution is not a corner solution
        r1 = NLsolve.nlsolve((x, fvec) -> f!(x, fvec, a_p, y_p, y_c, yi_p, yi_c),
                     [r.zero[1], 0.], iterations = 200)
        if converged(r1)
          r = r1
        end
      else
        a0 = max(0.,r.zero[1] + 3*randn(1)[1])
        b0 = max(0.,r.zero[2] + 3*randn(1)[1])
        if a0 + b0 > p.R * a_p + y_p
          a0 = p.R * a_p + y_p
          b  = 0
        end
        r = NLsolve.nlsolve((x, fvec) -> f!(x, fvec, a_p, y_p, y_c, yi_p, yi_c), [a0, b0], iterations = 100)
      end
      i  += 1
    end

    # ensure that the solution are not negative; should not happen
    if (!converged(r)) || (any(x -> x == true, r.zero .< -1e9))
      warn("Did not converged, ($a_p, $y_p, $y_c, $yi_p, $yi_c, $ai_p)")
      _A3, _B3 = pol.A3[ai_p,yi_p,yi_c], pol.B3[ai_p,yi_p,yi_c]
    else
      _A3, _B3 = r.zero
    end
    return (_A3, _B3)
  end

  for (ai_p, a_p) in enumerate(p.a_grid)
    for (yi_p, y_p) in enumerate(p.Y[:,3])
      for (yi_c,y_c) in enumerate(p.Y[:,1])
        _A3, _B3 = solver(a_p, y_p, y_c, yi_p, yi_c, ai_p)
        pol.A3[ai_p,yi_p,yi_c] = _A3 >= 0. ? _A3 : 0.
        pol.B3[ai_p,yi_p,yi_c] = _B3 >= 0. ? _B3 : 0.
      end
    end
  end
  info(" Old3! ended")
end

function old4!(p::Param, U::Utility, pol::Policies)
  A2, A3, B3 = interp(p, pol.A2, "A2"), interp(p, pol.A3, "A3"), interp(p, pol.B3, "B3")

  function f!(x, fvec, a_p::Float64, y_p::Float64, a_c::Float64, yi_c::Int, yi_p::Int)
    if x[1] + x[2] > p.R + a_p + y_p
      fvec[1] = 1e9
      fvec[2] = 1e9
    else
      Euc = 0
      _a2  = A2[yi_c][a_c, x[1], x[2]]
      # compute the t = 3 expectation
      for yi_fp in 1:length(p.Y[:,3])
        for yi_fc in 1:length(p.Y[:,1])
          _b3  = B3[yi_fp,yi_fc][_a2 + x[1]]
          _a3  = A3[yi_fp,yi_fc][_a2 + x[1]]
          Euc += U.du(p.R * (_a2 + x[1]) + p.Y[yi_fp,3]-_b3-_a3) * p.П[yi_c, yi_fp] * p.П[yi_c, yi_fc]
        end
      end
      fvec[1] = U.du(p.R * a_p + y_p - x[2] - x[1]) - p.α * p.R * p.β * Euc
      fvec[2] = U.du(p.R * a_p + y_p - x[2] - x[1]) - p.α * U.du(p.R * a_c + p.Y[yi_c,2] + x[2] - _a2)
    end
  end

  function solver(a_p::Float64, y_p::Float64, a_c::Float64, yi_c::Int, yi_p::Int)
    r = NLsolve.mcpsolve((x, fvec) -> f!(x, fvec, a_p, y_p, a_c, yi_c, yi_p), [0., 0.], [Inf, Inf],
                 [0., 0.], reformulation = :smooth, autodiff = true, iterations = 150)
    if !converged(r)
      # use the same procedure as for the old3! function. in the first iteration, always converges in
      # with the first mcpsolve
      r = NLsolve.nlsolve((x, fvec) -> f!(x, fvec, a_p, y_p, a_c, yi_c, yi_p),
                   r.zero, iterations = 150)
      if !converged(r)
        warn("Did not converged, ($a_p, $y_p, $a_c, $yi_c, $yi_p)")
      end
    end
    return r.zero
  end

  for (ai_p, a_p) in enumerate(p.a_grid)
    for (ai_c, a_c) in enumerate(p.a_grid)
      for (yi_p, y_p) in enumerate(p.Y[:,4])
        for yi_c in 1:length(p.Y[:,2])
          _A4, _B4 = solver(a_p, y_p, a_c, yi_c, yi_p)
          # Both may be negative due to approximations
          pol.A4[ai_p, ai_c, yi_p, yi_c] = _A4 >= 0. ? _A4 : 0.
          pol.B4[ai_p, ai_c, yi_p, yi_c] = _B4 >= 0. ? _B4 : 0.
        end
      end
    end
  end
  info(" Old4! ended")
end

function ite(maxIter::Int, p::Param;
             tol::Float64=1e-4, isCheck::Bool=true)

  U   = Utility(p)
  pol = Policies(p)

  for iter in 1:maxIter
    old4!(p, U, pol)
    young1!(p, U, pol)
    old3!(p, U, pol)
    err = young2!(p, U, pol)

    if isCheck
      info("Error @ iteration $iter: ", err)
    end
    if err < tol
      info("Find solutions after $iter iterations")
      return pol
      break
    elseif iter == maxIter
      error("No solutions found after $iter iterations")
      return pol
    end

  end

end


end