# uncomment and run to obtain the plots presented in the thesis

# using PyPlot, JLD, QuantEcon
# include("src/GenModel/PFI.jl")
# p      = GenPFI.Param(size=20, ut = "CRRA", ɣ = 3., α = 0.9);
# t      = load("output/env/perfect_pol_corr.jld");
# t      = t["pol"];
# pol_al = t[1];

# draw income shocks for an individual. starting values determined by parental income, Xp
function sim_shocks(p::GenPFI.Param,Xp::Int,isStart::Bool)
  # H = number of time periods
  H      = size(p.Y)[2]
  Ys     = size(p.Y)[1]
  # build a distribution for each time period
  P_dist = [[DiscreteRV(vec(p.П[i,:,j])) for i in 1:Ys] for j in 1:H]

  if isStart
    X      = Array{Int}(Int(H/2))
    for t in 1:Int(H/2)
      if t == 1
        X[t] = QuantEcon.draw(P_dist[3][Xp])
      else
        X[t] = QuantEcon.draw(P_dist[4][X[t-1]])
      end
    end
  else
    X      = Array{Int}(H)
    for t in 1:(H)
      if t == 1
        X[t] = QuantEcon.draw(P_dist[1][Xp])
      else
        X[t] = QuantEcon.draw(P_dist[t][X[t-1]])
      end
    end
  end
  return X
end

function sim_plot(p::GenPFI.Param,pol::GenPFI.Policies, gen::Int, a0::Float64, y0::Int)
  # matrix of shocks
  H = size(p.Y)[2]
  S = Array{Int}(gen, H)
  for g in 1:gen
    if g == 1
      S[g,:] = [0;y0;sim_shocks(p,y0,true)]
    else
      S[g,:] = sim_shocks(p,S[g-1,2],false)
    end
  end

  A3, B3 = GenPFI.interp_lin(p, pol.A3, "A3"), GenPFI.interp_lin(p, pol.B3, "B3");
  A1, A2 = GenPFI.interp_lin(p, pol.A1, "A1"), GenPFI.interp_lin(p, pol.A2, "A2");
  B4     = GenPFI.interp_lin(p, pol.B4, "B4");

  function sim_dec!(f::Int, S::Array{Int}, sav::Array{Float64}, tra::Array{Float64})
    # always simulate at the family level
    # start with (1,3) family
    yi_p, yi_c  = S[f,3], S[f+1,1]
    _a2         = sav[f, 2]
    _a3         = A3[yi_p,yi_c][p.R * _a2]
    _b3         = B3[yi_p,yi_c][p.R * _a2]
    _a3         = _a3 < 0. ? 0. : _a3
    _b3         = _b3 < 0. ? 0. : _b3
    _a1         = A1[yi_c,yi_p][_a3, _b3]
    _a1         = _a1 < 0. ? 0. : _a1
    sav[f,3]    = _a3
    tra[f,1]    = _b3
    cons[f,3]   = p.R * _a2 + p.Y[yi_p,3] - _b3 - _a3
    sav[f+1,1]  = _a1
    cons[f+1,1] = p.Y[yi_c,1] + _b3 - _a1

    # then: (2,4) family
    yi_p, yi_c = S[f,4], S[f+1,2]
    _b4 = B4[yi_p,yi_c][p.R * _a3,p.R * _a1]
    _b4 = _b4 < 0. ? 0. : _b4
    _a2 = A2[yi_c][p.R * _a1,_b4]
    _a2 = _a2 < 0. ? 0. : _a2
    tra[f,2]   = _b4
    cons[f,4]  = p.R * _a3 + p.Y[yi_p,4] - _b4
    sav[f+1,2] = _a2
    cons[f+1,2]= p.R * _a1 + p.Y[yi_c,2] + _b4 - _a2
  end

  sav      = zeros(gen,H-1)
  tra      = zeros(gen-1,Int(H/2))
  cons     = zeros(gen, H)
  aw       = zeros(2 * gen + 1)
  sav[1,2] = a0
  aw[1]    = a0
  t        = 1
  for f in 1:(gen-1)
    sim_dec!(f, S, sav, tra)
    aw[t + 1] = sav[f,3] + sav[f+1,1]
    aw[t + 2] = sav[f+1,2]
    t += 2
  end

  return (S, sav, tra, cons, aw)
end

function plot_sim(p::GenPFI.Param,pol::GenPFI.Policies, gen::Int, a0::Float64, y0::Int)

  S, sav, tra, cons, aw = sim_plot(p, pol, gen, a0, y0)

  # plot the results
  T      = 2 * gen + 1
  fig    = figure("simulation",figsize=(10,10))
  colors = get_cmap("Reds",gen)
  # subplots_adjust(wspace=0.4, hspace=0.4)
  # ax1    = subplot2grid((2, 2), (0, 0), colspan=2)
  subplots_adjust(hspace=0.0)
  ax1    = fig[:add_subplot](3,1,1)
  ax1[:plot](0:1, sav[1,2:3], linestyle = ":", marker = "o", color = colors(1))
  ax1[:plot](0:(length(aw)-2), aw[1:(end-1)], linestyle = "--", color = "black",
            alpha = 0.6, label = "Family wealth")
  ub = 1
  for g in 2:gen
    ax1[:plot](ub:(ub+2), sav[g,1:3], linestyle = ":", marker = "o", color = colors(g))
    ub = ub + 2
  end
  xlabel("Time")
  ylabel("Savings")
  # title("Savings")
  legend()

  # ax2    = subplot2grid((2, 2), (1, 0), colspan=1)
  ax2    = fig[:add_subplot](3,1,2,sharex=ax1)
  ax2[:plot](1:2, tra[1,1:2], linestyle = ":", marker = "o", color = colors(1))
  ub = 2
  for g in 2:(gen-1)
    ax2[:plot]((ub+1):(ub+2), tra[g,1:2], linestyle = ":", marker = "o", color = colors(g))
    ub = ub + 2
  end
  xlabel("Time")
  ylabel("Transfer")
  #title("Transfer")

  # ax3    = subplot2grid((2, 2), (1, 1), colspan=1)
  ax3    = fig[:add_subplot](3,1,3,sharex=ax1)
  ax3[:plot](0:2, [p.Y[S[1,i],i] for i in 2:4], linestyle = ":", marker = "o", color = colors(1))
  ub = 2
  for g in 2:gen
    ax3[:plot]((ub-1):(ub+2), [p.Y[S[g,i],i] for i in 1:4], linestyle = ":", marker = "o", color = colors(g))
    ub = (ub+2)
  end
  xlabel("Time")
  ylabel("Income")
  #title("Income shocks")
end
