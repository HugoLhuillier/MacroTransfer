
A2     = GenPFI.interp(p, pol.A2, "A2");
A3, B3 = GenPFI.interp(p, pol.A3, "A3"), GenPFI.interp(p, pol.B3, "B3");

function getParam4(ai_p::Int, ai_c::Int, yi_p::Int, yi_c::Int)
  return (p.a_grid[ai_p], p.a_grid[ai_c], p.Y[yi_p,4], p.Y[yi_c,2], yi_p, yi_c)
end

function f4(x::Array{Float64}, a_p::Float64, a_c::Float64, y_p::Float64, y_c::Float64, yi_p::Int, yi_c::Int)
  if x[1] + x[2] > p.R * a_p + y_p
    f1, f2 = 1e9, 1e9
  else
    _a2 = A2[yi_c][a_c,x[1],x[2]]
    if _a2 < 1e-3
      _a2, _da_a, _da_b = 0., 0., 0.
    else
      _da_a, _da_b = gradient(A2[yi_c], a_c, x[1], x[2])[2:3]
    end

    Eu3 = 0.
    for yi_fp in 1:length(p.Y[:,3])
      for yi_fc in 1:length(p.Y[:,1])
        _a3  = A3[yi_fp,yi_fc][x[1] + _a2]
        _b3  = B3[yi_fp,yi_fc][x[1] + _a2]
        _a3  = _a3 < 0. ? 0. : _a3
        _b3  = _b3 < 0. ? 0. : _b3
        Eu3 += U.du(p.R * (_a2 + x[1]) + p.Y[yi_fp,3] - _a3 - _b3) * p.П[yi_c,yi_fp,3] * p.П[yi_c,yi_fc,1]
      end
    end

    f1 = U.du(p.R * a_p + y_p - x[1] - x[2]) + p.α * _da_a * U.du(p.R * a_c + y_c + x[2] - _a2) -
              p.δ * p.R * (1 + _da_a) * Eu3
    f2 = U.du(p.R * a_p + y_p - x[1] - x[2]) - p.α * (1 - _da_b) * U.du(p.R * a_c + y_c + x[2] - _a2) -
              p.δ * p.R * _da_b * Eu3
  end
  return (f1, f2)
end

function foc4(par)
  n   = 20
  xx  = linspace(0.,par[1] + par[3]-0.1, n);
  yy  = linspace(0.,par[1] + par[3]-0.1, n);
  Z1  = zeros(n, n);
  Z2  = zeros(n, n);
  S   = zeros(n, n);

  for i in 1:n
    for j in 1:n
      Z1[j,i], Z2[j,i] = f4([xx[i],yy[j]],par...)
      if xx[i] + yy[j] > p.R * par[1] + par[3]
        S[j,i] = 1.
      else
        S[j,i] = NaN
      end
    end
  end

  xgrid = repmat(xx',n,1);
  ygrid = repmat(yy,1,n);
  fig = PyPlot.figure("FOC_cp",figsize=(10,10))
  ax = fig[:add_subplot](1,1,1)
  ax[:contourf](xgrid, ygrid, S, colors = "black", alpha = 0.3)
  cp1 = ax[:contour](xgrid, ygrid, Z2, [0.], colors = "red", linewidth = 3, zorder = 2)
  # ax[:clabel](cp1, inline=1, fontsize=10)
  cp2 = ax[:contour](xgrid, ygrid, Z1, [0.], linestyle = ":", colors = "blue", linewidth = 1,  zorder = 1)
  # ax[:clabel](cp2, inline=1, fontsize=10)
  PyPlot.xlabel(L"$a'$")
  PyPlot.ylabel(L"b")
end
