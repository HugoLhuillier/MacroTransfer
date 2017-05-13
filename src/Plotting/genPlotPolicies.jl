module GenPlotPolicies

include("./../GenModel/PFI.jl")
using PyPlot
info("Computation of the solution")
p   = GenPFI.Param(size = 5)
pol = GenPFI.ite(40,p;tol=4e-3)

function plot_4(isSurface::Bool)
  B4 = GenPFI.interp(p, pol.B4, "B4");
  n  = 50
  xx = linspace(0.,2*maximum(p.Y[:,4]), n);
  yy = linspace(0.,2*maximum(p.Y[:,2]), n);

  Zb1 = zeros(n, n);
  Zb2 = zeros(n, n);
  for i in 1:n
    for j in 1:n
      Zb1[j,i] = B4[1,4][xx[i], yy[j]]
      Zb2[j,i] = B4[4,1][xx[i], yy[j]]
    end
  end

  xgrid = repmat(xx',n,1);
  ygrid = repmat(yy,1,n);

  # with PyPlot
  fig = figure("policies_4",figsize=(10,10))
  if isSurface
    ax = fig[:add_subplot](1,2,1, projection = "3d")
    ax[:plot_surface](xgrid, ygrid, Zb, rstride=1, cstride=1, linewidth=0.1, cmap=ColorMap("bwr"))
    ax[:view_init](elev=14., azim=153.)
  else
    ax = fig[:add_subplot](1,1,1, projection = "3d")
    ax[:plot_wireframe](xgrid, ygrid, Zb1, rstride=5, cstride=5, linewidth=1, alpha = 0.8,
                        color = "salmon", label = L"$y_1, \tilde{y}_4$")
    ax[:plot_wireframe](xgrid, ygrid, Zb2, rstride=5, cstride=5, linewidth=1, alpha = 0.8,
                        color = "red", label = L"$y_4, \tilde{y}_1$")
    ax[:view_init](elev=11., azim=-106.)
    xlabel("Old's wealth")
    ylabel("Young's wealth")
    title("Age 4 transfer")
    legend(loc="center right")
  end

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
      Z1[j,i] = A1[1,1][xx[i], yy[j]]
      Z2[j,i] = A1[1,4][xx[i], yy[j]]
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
    ax[:plot_wireframe](xgrid, ygrid, Z1, rstride=5, cstride=5, linewidth=1,
                        color = "salmon", alpha = 0.8, label = L"$\tilde{y}_1, y_1$")
    ax[:plot_wireframe](xgrid, ygrid, Z2, rstride=5, cstride=5, linewidth=1,
                        color = "red", alpha = 0.8, label = L"$\tilde{y}_1, y_4$")
  end
  ax[:view_init](elev=26  ., azim=113.)
  xlabel("Old's savings")
  ylabel("Transfer received")
  title("Age 1 savings")
  legend(bbox_to_anchor=(0.05,0.5), loc=1)
end

function plot_2(isSurface)
  n  = 50
  A2 = GenPFI.interp(p, pol.A2, "A2");
  xx = linspace(0.,2*maximum(p.Y[:,3]), n);
  yy = linspace(0.,1/2*maximum(p.Y[:,3]), n);
  Z1  = zeros(n, n);
  Z2  = zeros(n, n);

  for i in 1:n
    for j in 1:n
      Z1[j,i] = A2[1][xx[i],yy[j]]
      Z2[j,i] = A2[4][xx[i],yy[j]]
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
    ax[:plot_wireframe](xgrid, ygrid, Z1, rstride=5, cstride=5, linewidth=1, alpha = 0.8,
                        color = "salmon", label = L"$\tilde{y}_1$")
    ax[:plot_wireframe](xgrid, ygrid, Z2, rstride=5, cstride=5, linewidth=1, alpha = 0.8,
                        color = "red", label = L"$\tilde{y}_4$")
  end
  ax[:view_init](elev=13., azim=-154.)
  xlabel("Young's savings")
  ylabel("Transfer received")
  title("Age 2 savings")
end

function plot_3_sav()
  A3 = GenPFI.interp(p, pol.A3, "A3")
  X  = linspace(0.,2*maximum(p.Y[:,3]),100)

  fig = figure("policies_3",figsize=(10,10))
  subplots_adjust(wspace=0.4, hspace=0.4)
  # ax1 = subplot2grid((1, 2), (0, 0), colspan=1)
  ax1 = fig[:add_subplot](2,2,1)
  ax1[:plot](X, A3[1,1][X], label = L"$\tilde{y}_1$", color = "pink")
  # ax1[:scatter](p.a_grid, A3[1,1][p.a_grid], label = "On grid", color = "black")
  ax1[:plot](X, A3[1,2][X], label = L"$\tilde{y}_2$", color = "salmon")
  ax1[:plot](X, A3[1,3][X], label = L"$\tilde{y}_3$", color = "indianred")
  ax1[:plot](X, A3[1,4][X], label = L"$\tilde{y}_4$", color = "darkred")
  xlabel("Wealth")
  ylabel(L"$\mathcal{A}_3$")
  title(L"$y_1$")
  legend(bbox_to_anchor=(-0.1, 0))

  ax2 = fig[:add_subplot](2,2,2)
  ax2[:plot](X, A3[2,1][X], color = "pink")
  # ax2[:scatter](p.a_grid, A3[2,1][p.a_grid], label = "On grid", color = "black")
  ax2[:plot](X, A3[2,2][X], color = "salmon")
  ax2[:plot](X, A3[2,3][X], color = "indianred")
  ax2[:plot](X, A3[2,4][X], color = "darkred")
  xlabel("Wealth")
  ylabel(L"$\mathcal{A}_3$")
  title(L"$y_2$")

  ax3 = fig[:add_subplot](2,2,3)
  ax3[:plot](X, A3[3,1][X], color = "pink")
  # ax1[:scatter](p.a_grid, A3[3,1][p.a_grid], label = "On grid", color = "black")
  ax3[:plot](X, A3[3,2][X], color = "salmon")
  ax3[:plot](X, A3[3,3][X], color = "indianred")
  ax3[:plot](X, A3[3,4][X], color = "darkred")
  xlabel("Wealth")
  ylabel(L"$\mathcal{A}_3$")
  title(L"$y_3$")

  ax4 = fig[:add_subplot](2,2,4)
  ax4[:plot](X, A3[4,1][X], color = "pink")
  # ax1[:scatter](p.a_grid, A3[4,1][p.a_grid], label = "On grid", color = "black")
  ax4[:plot](X, A3[4,2][X], color = "salmon")
  ax4[:plot](X, A3[4,3][X], color = "indianred")
  ax4[:plot](X, A3[4,4][X], color = "darkred")
  xlabel("Wealth")
  ylabel(L"$\mathcal{A}_3$")
  title(L"$y_4$")

end

function plot_3_trans()
  B3 = GenPFI.interp(p, pol.B3, "B3")
  X  = linspace(0.,2*maximum(p.Y[:,3]),100)

  fig = figure("policies_3",figsize=(10,10))
  subplots_adjust(wspace=0.4, hspace=0.4)
  ax1 = fig[:add_subplot](2,2,1)
  ax1[:plot](X, B3[1,1][X], label = L"$\tilde{y}_1$", color = "pink")
  ax1[:plot](X, B3[1,2][X], label = L"$\tilde{y}_2$", color = "salmon")
  ax1[:plot](X, B3[1,3][X], label = L"$\tilde{y}_3$", color = "indianred")
  ax1[:plot](X, B3[1,4][X], label = L"$\tilde{y}_4$", color = "darkred")
  xlabel("Wealth")
  ylabel(L"$\mathcal{B}_3$")
  title(L"$y_1$")
  legend(bbox_to_anchor=(-0.1, 0))

  ax2 = fig[:add_subplot](2,2,2)
  ax2[:plot](X, B3[2,1][X], color = "pink")
  ax2[:plot](X, B3[2,2][X], color = "salmon")
  ax2[:plot](X, B3[2,3][X], color = "indianred")
  ax2[:plot](X, B3[2,4][X], color = "darkred")
  xlabel("Wealth")
  ylabel(L"$\mathcal{B}_3$")
  title(L"$y_2$")

  ax3 = fig[:add_subplot](2,2,3)
  ax3[:plot](X, B3[3,1][X], color = "pink")
  ax3[:plot](X, B3[3,2][X], color = "salmon")
  ax3[:plot](X, B3[3,3][X], color = "indianred")
  ax3[:plot](X, B3[3,4][X], color = "darkred")
  xlabel("Wealth")
  ylabel(L"$\mathcal{B}_3$")
  title(L"$y_3$")

  ax4 = fig[:add_subplot](2,2,4)
  ax4[:plot](X, B3[4,1][X], color = "pink")
  # ax1[:scatter](p.a_grid, A3[4,1][p.a_grid], label = "On grid", color = "black")
  ax4[:plot](X, B3[4,2][X], color = "salmon")
  ax4[:plot](X, B3[4,3][X], color = "indianred")
  ax4[:plot](X, B3[4,4][X], color = "darkred")
  xlabel("Wealth")
  ylabel(L"$\mathcal{B}_3$")
  title(L"$y_4$")
end

end
