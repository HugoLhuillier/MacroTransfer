module GenPlotPolicies

include("./../GenModel/PFI.jl")
using PyPlot
info("Computation of the solution")
p   = GenPFI.Param(size=5)
pol = GenPFI.ite(40,p;tol=1e-6)

function plot_4(isSurface::Bool)
  B4 = GenPFI.interp(p, pol.B4, "B4");
  n  = 50
  xx = linspace(0.,2*maximum(p.Y[:,4]), n);
  yy = linspace(0.,2*maximum(p.Y[:,2]), n);

  Zb = zeros(n, n);
  for i in 1:n
    for j in 1:n
      Zb[j,i] = B4[1,1][xx[i], yy[j]]
    end
  end

  xgrid = repmat(xx',n,1);
  ygrid = repmat(yy,1,n);

  # with PyPlot
  fig = figure("policies_4",figsize=(10,10))
  ax = fig[:add_subplot](1,1,1, projection = "3d")
  if isSurface
    ax[:plot_surface](xgrid, ygrid, Zb, rstride=1, cstride=1, linewidth=0.1, cmap=ColorMap("bwr"))
  else
    ax[:plot_wireframe](xgrid, ygrid, Zb, rstride=5, cstride=5, linewidth=1, color = "red")
  end

  ax[:view_init](elev=14., azim=153.)
  xlabel("Old's wealth")
  ylabel("Young's wealth")
  title("Age 4 transfer")
end

function plot_1(isSurface)
  n  = 50
  A1 = GenPFI.interp(p, pol.A1, "A1");
  xx = linspace(0.,2*maximum(p.Y[:,3]), n);
  yy = linspace(0.,1/2*maximum(p.Y[:,3]), n);
  Z1  = zeros(n, n);
  Z2  = zeros(n, n);

  for i in 1:n
    for j in 1:n
      Z1[j,i] = A1[1,4][xx[i], yy[j]]
      Z2[j,i] = A1[4,1][xx[i], yy[j]]
    end
  end

  xgrid = repmat(xx',n,1);
  ygrid = repmat(yy,1,n);

  # with PyPlot
  fig = figure("policies_1",figsize=(10,10))
  ax = fig[:add_subplot](1,1,1, projection = "3d")
  if isSurface
    ax[:plot_surface](xgrid, ygrid, Z1, rstride=5, cstride=5, linewidth=1, cmap=ColorMap("Reds"))
  else
    ax[:plot_wireframe](xgrid, ygrid, Z1, rstride=5, cstride=5, linewidth=1, color = "red")
  end
  ax[:view_init](elev=18., azim=-119.)
  xlabel("Old's savings")
  ylabel("Transfer received")
  title("Age 1 savings")
end

function plot_2(isSurface)
  n  = 50
  A2 = GenPFI.interp(p, pol.A2, "A2");
  xx = linspace(0.,2*maximum(p.Y[:,3]), n);
  yy = linspace(0.,1/2*maximum(p.Y[:,3]), n);
  Z1  = zeros(n, n);

  for i in 1:n
    for j in 1:n
      Z1[j,i] = A2[1][xx[i],yy[j]]
    end
  end

  xgrid = repmat(xx',n,1);
  ygrid = repmat(yy,1,n);

  # with PyPlot
  fig = figure("policies_2",figsize=(10,10))
  ax = fig[:add_subplot](1,1,1, projection = "3d")
  if isSurface
    ax[:plot_surface](xgrid, ygrid, Z1, rstride=5, cstride=5, linewidth=1, cmap=ColorMap("Reds"))
  else
    ax[:plot_wireframe](xgrid, ygrid, Z1, rstride=5, cstride=5, linewidth=1, color = "red")
  end
  ax[:view_init](elev=13., azim=-154.)
  xlabel("Young's savings")
  ylabel("Transfer received")
  title("Age 2 savings")
end

function plot_3()
  A3, B3  = GenPFI.interp(p, pol.A3, "A3"), GenPFI.interp(p, pol.B3, "B3")
  X       = linspace(0.,2*maximum(p.Y[:,3]),100)

  fig = figure("policies_3",figsize=(10,10))
  ax = fig[:add_subplot](1,2,1)
  ax[:plot](X, A3[2,2][X], label = L"$y_2, \tilde{y}_2$", color = "pink")
  ax[:plot](X, A3[2,4][X], label = L"$y_2, \tilde{y}_4$", color = "salmon")
  ax[:plot](X, A3[4,2][X], label = L"$y_4, \tilde{y}_2$", color = "indianred")
  ax[:plot](X, A3[4,4][X], label = L"$y_4, \tilde{y}_4$", color = "darkred")
  xlabel("Wealth")
  ylabel(L"$\mathcal{A}_3$")
  title("Age 3 savings")
  legend()

  ax = fig[:add_subplot](1,2,2)
  ax[:plot](X, B3[2,2][X], label = L"$y_2, \tilde{y}_2$", color = "pink")
  ax[:plot](X, B3[2,4][X], label = L"$y_2, \tilde{y}_4$", color = "salmon")
  ax[:plot](X, B3[4,2][X], label = L"$y_4, \tilde{y}_2$", color = "indianred")
  ax[:plot](X, B3[4,4][X], label = L"$y_4, \tilde{y}_4$", color = "darkred")
  xlabel("Wealth")
  ylabel(L"$\mathcal{B}_3$")
  title("Age 3 transfer")
  legend()
end

end
