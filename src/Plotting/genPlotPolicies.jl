# uncomment and run to obtain the plots presented in the thesis

# using PyPlot, JLD
# include("src/GenModel/PFI.jl")
# p      = GenPFI.Param(size=20, ut = "CRRA", ɣ = 3., α = 0.9);
# t      = load("output/env/perfect_pol_corr.jld");
# t      = t["pol"];
# pol_al = t[1];

function plot_4(pol::GenPFI.Policies,p::GenPFI.Param;isSurface::Bool=false)
  B4 = GenPFI.interp_lin(p, pol.B4, "B4");
  n  = 50
  xx = linspace(0.,2*maximum(p.Y[:,4]), n);
  yy = linspace(0.,5*maximum(p.Y[:,2]), n);

  Zb1 = zeros(n, n);
  Zb2 = zeros(n, n);
  for i in 1:n
    for j in 1:n
      Zb1[j,i] = B4[1,2][xx[i], yy[j]]
      # Zb1[j,i] = B4[2,2][xx[i], yy[j]]
      Zb2[j,i] = B4[2,1][xx[i], yy[j]]
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
                        color = "salmon", label = L"$y_1, \tilde{y}_2$")
    ax[:plot_wireframe](xgrid, ygrid, Zb2, rstride=5, cstride=5, linewidth=1, alpha = 0.8,
                        color = "darkred", label = L"$y_2, \tilde{y}_1$")
    ax[:view_init](elev=11., azim=-106.)
    xlabel("Old's wealth")
    ylabel("Young's wealth")
    title("Age 4 transfer")
    legend(bbox_to_anchor=(0.05,0.5), loc=1)
  end

end

function plot_1(pol::GenPFI.Policies,p::GenPFI.Param;isSurface::Bool=false)
  n  = 50
  A1 = GenPFI.interp_lin(p, pol.A1, "A1");
  xx = linspace(0.,2*maximum(p.Y[:,3]), n);
  yy = linspace(0.,2*maximum(p.Y[:,3]), n);
  Z1  = zeros(n, n);
  Z2  = zeros(n, n);
  C1  = zeros(n,n);
  C2  = zeros(n,n);

  for i in 1:n
    for j in 1:n
      Z1[j,i] = A1[1,1][xx[i], yy[j]]
      # Z1[j,i] = A1[2,1][xx[i], yy[j]]
      Z2[j,i] = A1[2,2][xx[i], yy[j]]
      # Z2[j,i] = A1[1,2][xx[i], yy[j]]
      # C1[j,i] = p.Y[1,2] + yy[j] - A1[1,2][xx[i],yy[j]]
      # C2[j,i] = p.Y[2,2] + yy[j] - A1[1,2][xx[i],yy[j]]
    end
  end

  xgrid = repmat(xx',n,1);
  ygrid = repmat(yy,1,n);

  # with PyPlot
  fig = figure("policies_1",figsize=(10,10))
  # ax  = subplot2grid((2, 2), (0, 0), rowspan=2)
  ax = fig[:add_subplot](1,1,1, projection = "3d")
  if isSurface
    ax[:plot_surface](xgrid, ygrid, Z1, rstride=5, cstride=5, linewidth=1, cmap=ColorMap("Reds"))
  else
    ax[:plot_wireframe](xgrid, ygrid, Z1, rstride=5, cstride=5, linewidth=1,
                        color = "salmon", alpha = 0.8, label = L"$\tilde{y}_1, y_1$")
    ax[:plot_wireframe](xgrid, ygrid, Z2, rstride=5, cstride=5, linewidth=1,
                        color = "darkred", alpha = 0.8, label = L"$\tilde{y}_2, y_2$")
  end
  ax[:view_init](elev=18, azim=-19)
  xlabel("Old's savings")
  ylabel("Transfer received")
  title("Age 1 savings")
  legend(bbox_to_anchor=(0.05,0.5), loc=1)

end

function plot_2(pol::GenPFI.Policies,p::GenPFI.Param;isSurface::Bool=false)
  n  = 50
  A2 = GenPFI.interp_lin(p, pol.A2, "A2");
  xx = linspace(0.,2*maximum(p.Y[:,3]), n);
  yy = linspace(0.,2*maximum(p.Y[:,3]), n);
  Z1  = zeros(n, n);
  Z2  = zeros(n, n);
  C1  = zeros(n,n);

  for i in 1:n
    for j in 1:n
      Z1[j,i] = A2[1][xx[i],yy[j]]
      Z2[j,i] = A2[2][xx[i],yy[j]]

      # C1[j,i] = p.R * xx[i] + p.Y[1,2] + yy[j] - A2[1][xx[i],yy[j]]
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
                        color = "darkred", label = L"$\tilde{y}_2$")
  end
  ax[:view_init](elev=14., azim=134.)
  xlabel("Young's wealth")
  ylabel("Transfer received")
  title("Age 2 savings")
  legend(bbox_to_anchor=(0.15,0.5), loc=1)
end

function plot_3(pol::GenPFI.Policies,p::GenPFI.Param)
  A3 = GenPFI.interp_lin(p, pol.A3, "A3")
  B3 = GenPFI.interp_lin(p, pol.B3, "B3")
  # A3 = GenPFI.interp(p, pol.A3, "A3")
  # B3 = GenPFI.interp(p, pol.B3, "B3")
  Xs = linspace(0.,p.a_grid[end],100)

  fig = figure("policies_3",figsize=(10,10))
  subplots_adjust(wspace=0.4, hspace=0.4)
  # subplots_adjust(wspace=0.4, hspace=0.4)
  ax1 = fig[:add_subplot](2,2,1)
  ax1[:plot](Xs, A3[1,1][Xs], label = L"$\tilde{y}_1$", color = "salmon")
  # ax1[:scatter](p.a_grid, [pol.A3[i,1,1] for i in 1:length(p.a_grid)], color = "black")
  ax1[:plot](Xs, A3[1,2][Xs], label = L"$\tilde{y}_2$", color = "darkred")
  # ax1[:scatter](p.a_grid, [pol.A3[i,1,2] for i in 1:length(p.a_grid)], color = "black")
  xlabel("Wealth")
  ylabel(L"$\mathcal{A}_3$")
  title(L"Savings, $y_1$")
  legend()

  ax2 = fig[:add_subplot](2,2,2)
  ax2[:plot](Xs, A3[2,1][Xs], label = L"$\tilde{y}_1$", color = "salmon")
  # ax1[:scatter](p.a_grid, [pol.A3[i,1,1] for i in 1:length(p.a_grid)], color = "black")
  ax2[:plot](Xs, A3[2,2][Xs], label = L"$\tilde{y}_2$", color = "darkred")
  # ax1[:scatter](p.a_grid, [pol.A3[i,1,2] for i in 1:length(p.a_grid)], color = "black")
  xlabel("Wealth")
  ylabel(L"$\mathcal{A}_3$")
  title(L"Savings, $y_2$")
  legend()

  ax3 = fig[:add_subplot](2,2,3)
  ax3[:plot](Xs, B3[1,1][Xs], label = L"$\tilde{y}_1$", color = "salmon")
  # ax3[:scatter](p.a_grid, [pol.B3[i,1,1] for i in 1:length(p.a_grid)], color = "black")
  ax3[:plot](Xs, B3[1,2][Xs], label = L"$\tilde{y}_2$", color = "darkred")
  # ax3[:scatter](p.a_grid, [pol.B3[i,1,2] for i in 1:length(p.a_grid)], color = "black")
  xlabel("Wealth")
  ylabel(L"$\mathcal{B}_3$")
  title(L"Transfer, $y_1$")
  legend()

  ax3 = fig[:add_subplot](2,2,4)
  ax3[:plot](Xs, B3[2,1][Xs], label = L"$\tilde{y}_1$", color = "salmon")
  # ax3[:scatter](p.a_grid, [pol.B3[i,1,1] for i in 1:length(p.a_grid)], color = "black")
  ax3[:plot](Xs, B3[2,2][Xs], label = L"$\tilde{y}_2$", color = "darkred")
  # ax3[:scatter](p.a_grid, [pol.B3[i,1,2] for i in 1:length(p.a_grid)], color = "black")
  xlabel("Wealth")
  ylabel(L"$\mathcal{B}_3$")
  title(L"Transfer, $y_2$")
  legend()
end
