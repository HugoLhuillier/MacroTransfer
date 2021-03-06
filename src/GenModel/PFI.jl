module GenPFI

using Interpolations, ForwardDiff, NLsolve, Roots
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
  B4     = interp(p, pol.B4, "B4")
  B3, A3 = interp(p, pol.B3, "B3"), interp(p, pol.A3, "A3")
  A2     = interp(p, pol.A2, "A2")

  # FOC, to be solved for
  function f!(x, fvec, a_fp::Float64, b::Float64, y_c::Float64, yi_c::Int, yi_p::Int)
    # ensures that one cannot consume more than what one's have
    # less costly than putting an upper bound on the solver
    if (x[1] > y_c + b) || isequal(x[1],NaN)
      fvec[1] = 1e9
    else
      Eu2 = 0

      # compute the t = 2 expectation
      for yi_fc in 1:length(p.Y[:,2])
        for yi_fp in 1:length(p.Y[:,4])
          _b4  = B4[yi_fp, yi_fc][a_fp, x[1]]
          if _b4 > 1e-3
            # NOTE: gradient() do not access type ForwardDiff.Dual
            if typeof(x[1]) <: ForwardDiff.Dual{1,Float64}
              _db  = gradient(B4[yi_fp,yi_fc], a_fp, x[1].value)[2]
            else
              _db  = gradient(B4[yi_fp,yi_fc], a_fp, x[1])[2]
            end
          else
            _b4, _db = 0., 0.
          end
          _a2  = A2[yi_fc][x[1], _b4]
          _a2  = _a2 < 1e-3 ? 0. : _a2

          Eu2 += p.П[yi_c,yi_fc,2] * p.П[yi_p, yi_fp,4] *
                U.du(p.R * x[1] + p.Y[yi_fc, 2] + _b4 - _a2) * (p.R + _db)
        end
      end

      fvec[1] = U.du(y_c + b - x[1]) - p.β * Eu2
    end
  end

  function f(x, a_fp::Float64, b::Float64, y_c::Float64, yi_c::Int, yi_p::Int)
    Eu2 = 0
    Eu3 = 0

    # compute the t = 2 expectation
    for yi_fc in 1:length(p.Y[:,2])
      for yi_fp in 1:length(p.Y[:,4])
        _b4  = B4[yi_fp, yi_fc][a_fp, x]
        if _b4 > 1e-3
          _db  = gradient(B4[yi_fp,yi_fc], a_fp, x)[2]
        else
          _b4, _db = 0., 0.
        end
        _a2  = A2[yi_fc][x, _b4]
        _a2  = _a2 < 1e-3 ? 0. : _a2

        Eu2 += p.П[yi_c, yi_fc,2] * p.П[yi_p, yi_fp,4] *
              U.du(p.R * x + p.Y[yi_fc, 2] + _b4 - _a2) * (p.R + _db)
      end
    end

    return U.du(y_c + b - x) - p.β * Eu2
  end

  function solver(a_fp::Float64, b::Float64, y_c::Float64, yi_c::Int, yi_p::Int)
    # find the roots of the FOC
    # with non-monotonic function, mcpsolve will find that the boundary condition binds
    # only if start from 0. however, it also finds sometimes that it binds, even though
    # it does not (e.g. if f(0) > 0, f'(0) > 0 but f(.) = 0 exists due to non-monotonicity)
    # to circumvent this: (i) check that there is f(.) > 0 everywhere
    # (ii.a) if so, then corner binds
    # (ii.b) else, looking for an interior solution, with multiple starting values if necessary
    debug("New solver, ($a_fp, $b, $y_c, $yi_c, $yi_p)")

    X = linspace(0., b+y_c-1e-1,100)
    F = f.(X, a_fp,b,y_c,yi_c,yi_p)

    if any(x -> x == true, F .< 0)
      debug("Interior")
      r   = NLsolve.nlsolve((x, fvec) -> f!(x, fvec, a_fp, b, y_c, yi_c, yi_p),
                  [(b + y_c)/(1+p.R)], iterations = 150, autodiff = true)

      _a1 = r.zero[1]
      if !converged(r)
        debug("Not converged by middle")
        r   = NLsolve.nlsolve((x, fvec) -> f!(x, fvec, a_fp, b, y_c, yi_c, yi_p),
                      [0.], iterations = 100, autodiff = true)
        _a1 = r.zero[1]
        if !converged(r)
          debug("Not converged by 0")
          r = NLsolve.nlsolve((x, fvec) -> f!(x, fvec, a_fp, b, y_c, yi_c, yi_p),
                        [b + y_c - 1e-2], iterations = 100, autodiff = true)
          _a1 = r.zero[1]
          if !converged(r)
            # do brute force and use fzero from Roots package, more robust
            try
              _a1 = fzero(x -> f(x,a_fp,b,y_c,yi_c,yi_p),b+y_c-1e-1,[0.;b + y_c])
            catch
              warn("Did not converged, ($a_fp, $b, $y_c, $yi_c, $yi_p)")
              _a1 = 0.
            end
          end
        end
      end
    else
      debug("Corner")
      _a1 = 0.
    end

    return _a1
  end

  err = 0
  for (ai_fp, a_fp) in enumerate(p.a_grid)
    for (bi, b) in enumerate(p.b_grid)
      for (yi_c, y_c) in enumerate(p.Y[:,1])
        for (yi_p, y_p) in enumerate(p.Y[:,3])
          _A1 = solver(a_fp, b, y_c, yi_c, yi_p)
          _A1 = _A1 >= 1e-3 ? _A1 : 0.
          # _A1 may be negative due to the cubic interpolation
          err = abs(pol.A1[ai_fp,bi,yi_c,yi_p] - _A1) > err ?
                abs(pol.A1[ai_fp,bi,yi_c,yi_p] - _A1) : err
          pol.A1[ai_fp, bi, yi_c, yi_p] = _A1
        end
      end
    end
  end
  info("Young1! ended with error $err")

end

function young2!(p::Param, U::Utility, pol::Policies)
  A3, B3  = interp(p, pol.A3, "A3"), interp(p, pol.B3, "B3")

  function f!(x, fvec, a_c::Float64, b::Float64, y_c::Float64, yi_c::Int)
    # ensures that one cannot consume more than what one has
    if (x[1] > p.R * a_c + b + y_c) || (isequal(x[1],NaN))
      fvec[1] = 1e9
    else
      Eu3 = 0
      # compute the t = 2 expectation
      for yi_fp in 1:length(p.Y[:,3])
        for yi_fc in 1:length(p.Y[:,1])
          _b3 = B3[yi_fp,yi_fc][x[1]]
          _a3 = A3[yi_fp,yi_fc][x[1]]
          _b3 = _b3 < 1e-3 ? 0. : _b3
          _a3 = _a3 < 1e-3 ? 0. : _a3

          Eu3 += U.du(p.R * x[1] + p.Y[yi_fp,3] - _b3 - _a3) * p.П[yi_c,yi_fp,3] * p.П[yi_c,yi_fc,1]
        end
      end
      fvec[1] = U.du(p.R * a_c + b + y_c - x[1]) - p.R * p.β * Eu3
    end
  end

  function f(x, a_c::Float64, b::Float64, y_c::Float64, yi_c::Int)
    # ensures that one cannot consume more than what one has
    if x > p.R * a_c + b + y_c
      return 1e9
    else
      Eu3 = 0
      # compute the t = 2 expectation
      for yi_fp in 1:length(p.Y[:,3])
        for yi_fc in 1:length(p.Y[:,1])
          _b3 = B3[yi_fp,yi_fc][x]
          _a3 = A3[yi_fp,yi_fc][x]
          _b3 = _b3 < 1e-3 ? 0. : _b3
          _a3 = _a3 < 1e-3 ? 0. : _a3

          Eu3 += U.du(p.R * x + p.Y[yi_fp,3] - _b3 - _a3) * p.П[yi_c,yi_fp,3] * p.П[yi_c,yi_fc,1]
        end
      end
      return U.du(p.R * a_c + b + y_c - x) - p.R * p.β * Eu3
    end
  end

  function solver(a_c::Float64, b::Float64, y_c::Float64, yi_c::Int)
    # find the root of the FOC. use same method as young!1, see above
    X = linspace(0., p.R * a_c + b + y_c - 1e-1,100)
    F = f.(X, a_c,b,y_c,yi_c)
    if any(x -> x == true, F .< 0)
      debug("Interior")
      r   = NLsolve.nlsolve((x, fvec) -> f!(x, fvec, a_c, b, y_c, yi_c),
                   [(p.R * a_c + b + y_c) / (1 + p.R)], iterations = 150, autodiff = true)
      _a2 = r.zero[1]
      if !converged(r)
        debug("Not converged by middle")
        r   = NLsolve.nlsolve((x, fvec) -> f!(x, fvec, a_c, b, y_c, yi_c),
                     [0.], iterations = 150, autodiff = true)
        _a2 = r.zero[1]
        if !converged(r)
          debug("Not converged by 0")
          r   = NLsolve.nlsolve((x, fvec) -> f!(x, fvec, a_c, b, y_c, yi_c),
                        [p.R * a_c + b + y_c - 1e-2], iterations = 150, autodiff = true)
          _a2 = r.zero[1]
          if !converged(r)
            try
              _a2 = fzero(x -> f(x,a_c,b,y_c,yi_c),p.R * a_c +b+y_c-1e-1,[0.;p.R * a_c +b + y_c])
            catch
              warn("Did not converged, ($b, $y_c, $yi_c)")
              _a2 = 0.
            end
          end
        end
      end
    else
      debug("Corner")
      _a2 = 0.
    end
    return _a2
  end

  err = 0
  for (ai_c, a_c) in enumerate(p.a_grid)
      for (bi, b) in enumerate(p.b_grid)
        for (yi_c, y_c) in enumerate(p.Y[:,2])
          _A2 = solver(a_c,b,y_c,yi_c)
          _A2 = _A2 >= 1e-3 ? _A2 : 0.
          # update the distance between the two policy functions
          err = abs(pol.A2[ai_c,bi,yi_c] - _A2) > err ?
                abs(pol.A2[ai_c,bi,yi_c] - _A2) : err
          pol.A2[ai_c,bi,yi_c] = _A2
        end
    end
  end
  info(" Young2! ended")

  return err
end

function old3!(p::Param, U::Utility, pol::Policies, init::Bool)
  # markov strategies to be interpolated
  B4      = interp(p, pol.B4, "B4")
  A1      = interp(p, pol.A1, "A1")
  A2      = interp(p, pol.A2, "A2")

  function f!(x, fvec, a_p::Float64, y_p::Float64, y_c::Float64, yi_p::Int, yi_c::Int)
    if x[1] + x[2] > p.R * a_p + y_p || any(x -> x == true, isequal.(x,NaN))
      fvec[1] = 1e9
      fvec[2] = 1e9
    else
      _a1 = A1[yi_c, yi_p][x[1], x[2]]
      if _a1 < 1e-3
        _a1          = 0.
        _da_a, _da_b = 0., 0.
      else
        if typeof(x[1]) <: ForwardDiff.Dual{2,Float64}
          _da_a, _da_b = gradient(A1[yi_c, yi_p], x[1].value, x[2].value)
        else
          _da_a, _da_b = gradient(A1[yi_c, yi_p], x[1], x[2])
        end
      end

      # compute the t = 4 expectations
      Eu4, Eu2 = 0., 0.
      for yi_fp in 1:length(p.Y[:,4])
        for yi_fc in 1:length(p.Y[:,2])
          _b4  = B4[yi_fp,yi_fc][x[1],_a1]
          _b4  = _b4 < 1e-3 ? 0. : _b4
          _a2  = A2[yi_fc][_a1,_b4]
          _a2  = _a2 < 1e-3 ? 0. : _a2

          Eu4 += U.du(p.R * x[1] + p.Y[yi_fp,4] - _b4) * p.П[yi_p,yi_fp,4] * p.П[yi_c,yi_fc,2]
          Eu2 += U.du(p.R * _a1 + p.Y[yi_fc,2] + _b4 - _a2) * p.П[yi_p,yi_fp,4] * p.П[yi_c,yi_fc,2]
        end
      end
      fvec[1] = U.du(p.R * a_p + y_p - x[2] - x[1]) + p.α * _da_a * U.du(y_c + x[2] - _a1) -
              p.R * p.β * (Eu4 + p.α * _da_a * Eu2)
      fvec[2] = U.du(p.R * a_p + y_p - x[2] - x[1]) - p.α * (1 - _da_b) * U.du(y_c + x[2] - _a1) -
              p.α * p.R * p.β * _da_b * Eu2
    end
  end

  function f_spe!(x, fvec, a_p::Float64, y_p::Float64, y_c::Float64, yi_p::Int, yi_c::Int)
    if x[1] > p.R * a_p + y_p || any(x -> x == true, isequal.(x,NaN))
      fvec[1] = 1e9
    else
      _a1 = A1[yi_c, yi_p][x[1], 0.]
      if _a1 < 1e-3
        _a1          = 0.
        _da_a, _da_b = 0., 0.
      else
        if typeof(x[1]) <: ForwardDiff.Dual{2,Float64}
          _da_a, _da_b = gradient(A1[yi_c, yi_p], x[1].value, 0.)
        else
          _da_a, _da_b = gradient(A1[yi_c, yi_p], x[1], 0.)
        end
      end

      # compute the t = 4 expectations
      Eu4, Eu2 = 0., 0.
      for yi_fp in 1:length(p.Y[:,4])
        for yi_fc in 1:length(p.Y[:,2])
          _b4  = B4[yi_fp,yi_fc][x[1],_a1]
          _b4  = _b4 < 1e-3 ? 0. : _b4
          _a2  = A2[yi_fc][_a1,_b4]
          _a2  = _a2 < 1e-3 ? 0. : _a2

          Eu4 += U.du(p.R * x[1] + p.Y[yi_fp,4] - _b4) * p.П[yi_p,yi_fp,4] * p.П[yi_c,yi_fc,2]
          Eu2 += U.du(p.R * _a1 + p.Y[yi_fc,2] + _b4 - _a2) * p.П[yi_p,yi_fp,4] * p.П[yi_c,yi_fc,2]
        end
      end
      fvec[1] = U.du(p.R * a_p + y_p - x[1]) + p.α * _da_a * U.du(y_c - _a1) -
              p.R * p.β * (Eu4 + p.α * _da_a * Eu2)
    end
  end

  function f(x, a_p::Float64, y_p::Float64, y_c::Float64, yi_p::Int, yi_c::Int)
    if x[1] + x[2] > p.R * a_p + y_p
      f1, f2 = 1e9, 1e9
    else
      _a1            = A1[yi_c, yi_p][x[1], x[2]]
      if _a1 < 1e-3
        _a1          = 0.
        _da_a, _da_b = 0., 0.
      else
        _da_a, _da_b = gradient(A1[yi_c, yi_p], x[1], x[2])
      end
      # compute the t = 4 expectations
      Eu4, Eu2 = 0., 0.
      for yi_fp in 1:length(p.Y[:,4])
        for yi_fc in 1:length(p.Y[:,2])
          _b4  = B4[yi_fp,yi_fc][x[1],_a1]
          _b4  = _b4 < 1e-3 ? 0. : _b4
          _a2  = A2[yi_fc][_a1,_b4]
          _a2  = _a2 < 1e-3 ? 0. : _a2

          Eu4 += U.du(p.R * x[1] + p.Y[yi_fp,4] - _b4) * p.П[yi_p,yi_fp,4] * p.П[yi_c,yi_fc,2]
          Eu2 += U.du(p.R * _a1 + p.Y[yi_fc,2] + _b4 - _a2) * p.П[yi_p,yi_fp,4] * p.П[yi_c,yi_fc,2]
        end
      end
      f1 = U.du(p.R * a_p + y_p - x[2] - x[1]) + p.α * _da_a * U.du(y_c + x[2] - _a1) -
              p.R * p.β * (Eu4 + p.α * _da_a * Eu2)
      f2 = U.du(p.R * a_p + y_p - x[2] - x[1]) - p.α * (1 - _da_b) * U.du(y_c + x[2] - _a1) -
              p.α * p.R * p.β * _da_b * Eu2
    end
    return (f1,f2)
  end

  function solver(a0::Float64, b0::Float64,
                  a_p::Float64, y_p::Float64, y_c::Float64, yi_p::Int, yi_c::Int, ai_p::Int)
    # find the roots of the system of equations, for given starting values
    r = NLsolve.mcpsolve((x, fvec) -> f!(x, fvec, a_p, y_p, y_c, yi_p, yi_c), [0., 0.], [Inf, Inf],
                 [a0, b0], reformulation = :smooth, iterations = 150, autodiff = true)
    # in most of the cases, the solution does exist, but mcpsolve did not find it => multi-start approach
    i = 1
    while (!converged(r)) && (i <= 10)
      # NOTE: to find the roots, we use instead nlsolve, as it works better than mcpsolve
      # also, note that we can be sure that the solution is not a corner solution
      a0 = max(0.,r.zero[1] + 2*randn(1)[1])
      b0 = max(0.,r.zero[2] + 2*randn(1)[1])
      if a0 + b0 > p.R * a_p + y_p
        a0 = p.R * a_p + y_p
        b  = 0
      end
      r = NLsolve.nlsolve((x, fvec) -> f!(x, fvec, a_p, y_p, y_c, yi_p, yi_c),
                          [a0, b0], iterations = 100, autodiff = true)
      i  += 1
    end

    # ensure that the solution are not negative; should not happen
    if (!converged(r)) || (any(x -> x == true, r.zero .< -1e9))
      _a3, _b3 = Void, Void
    else
      _a3, _b3 = r.zero
    end
    return (_a3, _b3)
  end

  function pick_root(r::Array)
    _a1, _b1 = 0., 0.
    _a2, _b2 = 0., 0.
    for i in 1:size(r)[1]
      if r[i][1] > _a1
        _a1, _b1 = r[i]
      end
      if r[i][2] > _b2
        _a2, _b2 = r[i]
      end
    end
    if (_a1 < _a2) || (_b1 > _b2)
      warn("Did not work")
    end
    _a = _a2 + (1/2) * (_a1 - _a2)
    _b = _b1 + (1/2) * (_b2 - _b1)

    return (_a, _b)
  end

  function find_root(maxRoots::Int, a_p::Float64, y_p::Float64, y_c::Float64,
                    yi_p::Int, yi_c::Int, ai_p::Int)

    start    = [(pol.A3[ai_p,yi_p,yi_c], pol.B3[ai_p,yi_p,yi_c]); (0.,0.);
                (0.,p.R * a_p + y_p);(p.R * a_p + y_p,0.);
                ((p.R * a_p + y_p)/2, (p.R * a_p + y_p)/2)]

    for i in size(start)[1]:maxRoots
      a0 = y_p * rand(1)[1]
      b0 = p.α * y_p * rand(1)[1]
      if a0 + b0 > p.R * a_p + y_p
        a0 = p.R * a_p + y_p - 1e-1
        b0 = 0
      end
      push!(start,(a0,b0))
    end

    r        = Tuple[]
    for i in 1:maxRoots
      _a, _b = solver(start[i]..., a_p,y_p,y_c,yi_p,yi_c,ai_p)
      if all(x -> x == true, [_a,_b] .!= Void)
        _a, _b = round(_a,3), round(_b,3)
        if _a < 0.
          _a = 0.
        end
        if _b < 0.
          _b = 0.
        end
        if all(x -> x == false, [isequal(r[i],(_a,_b)) for i in 1:size(r)[1]])
          push!(r, (_a, _b))
          debug("New root: ($_a, $_b)")
        end
      end
    end

    if size(r)[1] == 0
      return (0., 0.)
    else
      return pick_root(r)
    end
  end

  function find_root_spe(maxRoots::Int, a_p::Float64, y_p::Float64, y_c::Float64,
                    yi_p::Int, yi_c::Int, ai_p::Int)

    r   = Tuple[]
    #1. look for the root with b = 0.; if f(a,0.) = 0, good
    n   = NLsolve.nlsolve((x, fvec) -> f_spe!(x, fvec, a_p,y_p,y_c,yi_p,yi_c),
                 [(p.R * a_p + y_p) / (1 + p.R)], iterations = 150, autodiff = true)
    if !converged(n)
      debug("Not converged by middle")
      n   = NLsolve.nlsolve((x, fvec) -> f!(x, fvec, a_c, b, y_c, yi_c),
                  [0.], iterations = 150, autodiff = true)
      if !converged(n)
        debug("Not converged by 0")
        n   = NLsolve.nlsolve((x, fvec) -> f!(x, fvec, a_c, b, y_c, yi_c),
                     [p.R * a_p + y_p - 1e-2], iterations = 150, autodiff = true)
      end
    end
    a0, b0 = n.zero[1], 0.
    f1, f2 = f([a0,b0], a_p, y_p, y_c, yi_p, yi_c)
    cor    = false
    if (abs(f1) < 1e-4) && (abs(f2) < 1e-4)
      # info("Found one sol for b = 0 & ($ai_p, $yi_p, $yi_c)")
      push!(r, (a0, b0))
      cor  = true
    end

    #2. Look for any other root
    start    = [(pol.A3[ai_p,yi_p,yi_c], pol.B3[ai_p,yi_p,yi_c]); (0.,0.);
                (0.,p.R * a_p + y_p);(p.R * a_p + y_p,0.);
                ((p.R * a_p + y_p)/2, (p.R * a_p + y_p)/2)]

    for i in size(start)[1]:maxRoots
      a0 = y_p * rand(1)[1]
      b0 = p.α * y_p * rand(1)[1]
      if a0 + b0 > p.R * a_p + y_p
        a0 = p.R * a_p + y_p - 1e-1
        b0 = 0
      end
      push!(start,(a0,b0))
    end

    for i in 1:maxRoots
      _a, _b = solver(start[i]..., a_p,y_p,y_c,yi_p,yi_c,ai_p)
      if all(x -> x == true, [_a,_b] .!= Void)
        _a, _b = round(_a,3), round(_b,3)
        if (all(x -> x == false, [isequal(r[i],(_a,_b)) for i in 1:size(r)[1]])) && (_b > 0.)
          push!(r, (_a, _b))
          debug("New root: ($_a, $_b)")
        end
      end
    end

    #3. and pick the one with the highest b (and highest a if 1 failed)
    function pick_root_spe(r::Array)
      _a1, _b1 = 0., 0.
      _a2, _b2 = 0., 0.
      for i in 1:size(r)[1]
        if r[i][2] == 0.
          _a1, _b1 = r[i]
        end
        if r[i][2] > _b2
          _a2, _b2 = r[i]
        end
      end
      if (_a1 < _a2) || (_b1 > _b2)
        warn("Did not work")
      end
      _a = _a2 + (1/2) * (_a1 - _a2)
      _b = _b1 + (1/2) * (_b2 - _b1)
      return (_a, _b)
    end

    if size(r)[1] == 0
      # warn("Here @ ($ai_p, $yi_p, $yi_c)")
      return (0., 0.)
    elseif !cor
      return pick_root(r)
    else
      return pick_root_spe(r)
    end
  end

  function gen_solver(a_p::Float64,y_p::Float64,y_c::Float64,yi_p::Int,yi_c::Int,
                      ai_p::Int, init::Bool)

    n   = 40
    xx  = linspace(0.,p.R * a_p + y_p, n);
    yy  = linspace(0.,p.R * a_p + y_p, n);
    F1  = zeros(n, n);
    F2  = zeros(n, n);

    for i in 1:n
      for j in 1:n
        F1[j,i], F2[j,i] = f([xx[i],yy[j]],a_p,y_p,y_c,yi_p,yi_c)
      end
    end

    if (all(x -> x == true, F1 .> 0)) && (all(x -> x == true, F2 .> 0))
      _a, _b = 0., 0.
      # info("Corner solution")
    else
      if !init
        # if (yi_p == 1) && (yi_c == 2)
          _a, _b = find_root_spe(100, a_p, y_p, y_c, yi_p, yi_c, ai_p)
          i      = 1
          if ai_p > 1
            # info("_a = $_a vs. $(pol.A3[ai_p-1,yi_p,yi_c]) && _b = $_b vs. $(pol.B3[ai_p-1,yi_p,yi_c])")
            a_b, b_b = pol.A3[ai_p-1,yi_p,yi_c], pol.B3[ai_p-1,yi_p,yi_c]
            while ((_a + 1e-1 < a_b) || (_b + 1e-1 < b_b))
              _a, _b = find_root_spe(50, a_p, y_p, y_c, yi_p, yi_c, ai_p)
              i      += 1
              if i == 10
                warn("Did not find the correct root @ ($ai_p, $yi_p, $yi_c)")
                break
              end
            end
          end
      else
        _a, _b = find_root(100, a_p, y_p, y_c, yi_p, yi_c, ai_p)
      end
    end

    return (_a, _b)
  end

  err = 0
  for (yi_p, y_p) in enumerate(p.Y[:,3])
    for (yi_c,y_c) in enumerate(p.Y[:,1])
      for (ai_p, a_p) in enumerate(p.a_grid)
        _A3, _B3 = gen_solver(a_p, y_p, y_c, yi_p, yi_c, ai_p, init)
        _A3 = _A3 >= 1e-3 ? _A3 : 0.
        _B3 = _B3 >= 1e-3 ? _B3 : 0.
        # println("@($ai_p, $yi_p, $yi_c) \n A: $(pol.A3[ai_p,yi_p,yi_c]) vs. $_A3 \n B: $(pol.B3[ai_p,yi_p,yi_c]) vs. $_B3")
        err = abs(pol.A3[ai_p,yi_p,yi_c] - _A3) > err ?
              abs(pol.A3[ai_p,yi_p,yi_c] - _A3) : err
        err = abs(pol.B3[ai_p,yi_p,yi_c] - _B3) > err ?
              abs(pol.B3[ai_p,yi_p,yi_c] - _B3) : err
        pol.A3[ai_p,yi_p,yi_c] = _A3
        pol.B3[ai_p,yi_p,yi_c] = _B3 
      end
    end
  end
  info(" Old3! ended with error $err")
end

function old4!(p::Param, U::Utility, pol::Policies)
  A2, A3, B3 = interp(p, pol.A2, "A2"), interp(p, pol.A3, "A3"), interp(p, pol.B3, "B3")

  function f!(x, fvec, a_p::Float64, y_p::Float64, a_c::Float64, yi_c::Int, yi_p::Int)
    if x[1] > p.R + a_p + y_p || isequal(x[1],NaN)
      fvec[1] = 1e9
    else
      Euc  = 0
      _a2  = A2[yi_c][a_c, x[1]]
      _a2  = _a2 <1e-3 ? 0. : _a2
      # compute the t = 3 expectation
      for yi_fp in 1:length(p.Y[:,3])
        for yi_fc in 1:length(p.Y[:,1])
          _b3  = B3[yi_fp,yi_fc][_a2]
          _a3  = A3[yi_fp,yi_fc][_a2]
          _b3  = _b3 < 1e-3 ? 0. : _b3
          _a3  = _a3 < 1e-3 ? 0. : _a3

          Euc += U.du(p.R * _a2 + p.Y[yi_fp,3] - _b3 - _a3) * p.П[yi_c, yi_fp,3] * p.П[yi_c, yi_fc,1]
        end
      end
      fvec[1] = U.du(p.R * a_p + y_p - x[1]) - p.α * U.du(p.R * a_c + p.Y[yi_c,2] + x[1] - _a2)
    end
  end

  function f(x, a_p::Float64, y_p::Float64, a_c::Float64, yi_c::Int, yi_p::Int)
    if x > p.R + a_p + y_p
      f1 = 1e9
    else
      Euc  = 0
      _a2  = A2[yi_c][a_c, x]
      _a2  = _a2 < 1e-3 ? 0. : _a2
      # compute the t = 3 expectation
      for yi_fp in 1:length(p.Y[:,3])
        for yi_fc in 1:length(p.Y[:,1])
          _b3  = B3[yi_fp,yi_fc][_a2]
          _a3  = A3[yi_fp,yi_fc][_a2]
          _b3  = _b3 < 1e-3 ? 0. : _b3
          _a3  = _a3 < 1e-3 ? 0. : _a3

          Euc += U.du(p.R * _a2 + p.Y[yi_fp,3]-_b3-_a3) * p.П[yi_c, yi_fp,3] * p.П[yi_c, yi_fc,1]
        end
      end
      f1 = U.du(p.R * a_p + y_p - x) - p.α * U.du(p.R * a_c + p.Y[yi_c,2] + x - _a2)
    end
    return f1
  end

  function solver(a_p::Float64, y_p::Float64, a_c::Float64, yi_c::Int, yi_p::Int)
    # NOTE: takes into account the multiplicity of equilibrium
    # check whether both choices are going to lie in the interior set
    B  = linspace(0.,p.R * a_p + y_p,100)
    F1 = f.(B,a_p,y_p,a_c,yi_c,yi_p)

    if (any(x -> x == true, F1 .< 0))
      # there exists an interior solution for both choices => set x[1] -> 0, and max on x[2]
      # _a = 1e-7
      r  = NLsolve.mcpsolve((x, fvec) -> f!(x, fvec, a_p, y_p, a_c, yi_c, yi_p), [0.], [Inf],
                   [p.R * a_p + y_p - 1], reformulation = :smooth, autodiff = true, iterations = 100)
      _b = r.zero[1]
      if !converged(r)
        r  = NLsolve.nlsolve((x, fvec) -> f!(x, fvec, a_p, y_p, a_c, yi_c, yi_p),
                     [r.zero[1]-1], iterations = 100, autodiff = true)
        _b = r.zero[1]
        if !converged(r)
          try
            _b = fzero(x -> f(x,a_p,y_p,a_c,yi_c,yi_p),p.R * a_p + y_p-1,[0.;p.R * a_p + y_p])
          catch
            warn("Did not converged, ($a_fp, $b, $y_c, $yi_c, $yi_p)")
            _b = 0.
          end
        end
      end
    else
      _b = 0.
    end
    return _b
  end

  err = 0
  for (ai_p, a_p) in enumerate(p.a_grid)
    for (ai_c, a_c) in enumerate(p.a_grid)
      for (yi_p, y_p) in enumerate(p.Y[:,4])
        for yi_c in 1:length(p.Y[:,2])
          _B4 = solver(a_p, y_p, a_c, yi_c, yi_p)
          _B4 = _B4 >= 1e-3 ? _B4 : 0.
          # Both may be negative due to approximations
          err = abs(pol.B4[ai_p,ai_c,yi_p,yi_c] - _B4) > err ?
                abs(pol.B4[ai_p,ai_c,yi_p,yi_c] - _B4) : err
          pol.B4[ai_p, ai_c, yi_p, yi_c] = _B4
        end
      end
    end
  end
  info(" Old4! ended with error $err")
end

function ite(maxIter::Int, p::Param;
             tol::Float64=1e-3, isCheck::Bool=true)

  U   = Utility(p)
  pol = Policies(p)
  err = 20.

  for iter in 1:maxIter
    old4!(p, U, pol)
    young1!(p, U, pol)
    if err > 1
      old3!(p, U, pol,true)
    else
      old3!(p, U, pol,false)
    end
    err = young2!(p, U, pol)

    if isCheck
      println("\t\tError @ iteration $iter: ", err)
    end
    if err < tol
      info("Find solutions after $iter iterations")
      return pol
      break
    elseif iter == maxIter
      warn("No solutions found after $iter iterations")
      return pol
    end

  end

end


end
