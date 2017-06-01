A4, B4 = GenPFI.interp(p, pol.A4, "A4"), GenPFI.interp(p, pol.B4, "B4");
A1, A2 = GenPFI.interp(p, pol.A1, "A1"), GenPFI.interp(p, pol.A2, "A2");

function getParam3(ai_p, yi_p, yi_c)
  return (p.a_grid[ai_p], p.Y[yi_p,3], p.Y[yi_c,1], yi_p, yi_c)
end

function f3(x::Array{Float64}, a_p::Float64, y_p::Float64, y_c::Float64, yi_p::Int, yi_c::Int)
  if x[1] + x[2] > p.R * a_p + y_p
    f1, f2 = 1e9, 1e9
  else
    _a1 = A1[yi_c,yi_p][x[1],x[2]]
    if _a1 < 0.
      _a1, _da_a, _da_b = 0., 0., 0.
    else
      _da_a, _da_b = gradient(A1[yi_c, yi_p], x[1], x[2])
    end

    Eu4, Eu2 = 0., 0.
    for yi_fp in 1:length(p.Y[:,4])
      for yi_fc in 1:length(p.Y[:,2])
        _a4 = A4[yi_fp,yi_fc][x[1], _a1]
        _b4 = B4[yi_fp,yi_fc][x[1], _a1]
        _a4 = _a4 < 0. ? 0. : _a4
        _b4 = _b4 < 0. ? 0. : _b4
        _a2 = A2[yi_fc][_a1, _a4, _b4]
        _a2 = _a2 < 0. ? 0. : _a2

        Eu4 += U.du(p.R * x[1] + p.Y[yi_fp,4] - _b4 - _a4) * p.П[yi_p,yi_fp,4] * p.П[yi_c,yi_fc,2]
        Eu2 += U.du(p.R * _a1 + p.Y[yi_fc,2] + _b4 - _a2) * p.П[yi_p,yi_fp,4] * p.П[yi_c,yi_fc,2]
      end
    end

    f1 = U.du(p.R * a_p + y_p - x[1] - x[2]) - p.α * _da_a * U.du(y_c + x[2] - _a1) -
              p.β * p.R * (Eu4 + p.α * _da_a * Eu2)
    f2 = U.du(p.R * a_p + y_p - x[1] - x[2]) - p.α * (1 - _da_b) * U.du(y_c + x[2] - _a1) -
              p.α * p.β * p.R * _da_b * Eu2
  end
  return (f1, f2)
end

function foc3(par)
  n   = 20
  xx  = linspace(0.,par[1] + par[2], n);
  yy  = linspace(0.,par[1] + par[2], n);
  Z1  = zeros(n, n);
  Z2  = zeros(n, n);
  S   = zeros(n, n);

  for i in 1:n
    for j in 1:n
      Z1[j,i], Z2[j,i] = f3([xx[i],yy[j]],par...)
      if xx[i] + yy[j] > p.R * par[1] + par[2]
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
