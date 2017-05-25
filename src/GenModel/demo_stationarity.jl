using StatsBase, QuantEcon
using DataFrames
include("src/GenModel/PFI.jl")
include("src/GenModel/stationary_dist.jl")
using Plots
using JLD

p      = GenPFI.Param(size=20, ut = "CRRA", ɣ = 3., α = 0.9);
# load the policy function
t      = load("output/env/perfect_pol_corr.jld");
t      = t["pol"];
pol_al = t[1];
# to directly used the stationary distributions of the thesis, run these
t      = load("output/env/fam_al.jld");
t      = t["fam"];
fam_al = t[1]
t      = load("output/env/Fam_noal.jld");
t      = t["fam"];
fam_noal = t[1];
# otherwise, can also do (here for the altruistic case)
fam_al = stationarity(p, pol_al, 500_000)
# simulate a two-period panel from the stationary distribution
data   = simulate_pan(2, fam_al, p, pol_al)
# .. and add manually some more variables
data[:trans_bin]   = Int.(map((x) -> x > 0 ? 1 : 0, data[:b]));
data[:dummy_age]   = Int.(map(x -> x == 3 ? 0 : 1, data[:g]));
data[:constrained] = Int.(map(x -> x < 0.1 ? 1 : 0, data[:a_c]));
data[:c_o] = (map((x,y,z,w) -> p.R * x + y - w - z, data[:a_p], data[:y_p], data[:a_fp], data[:b]));
data[:c_y] = (map((x,y,w,z) -> p.R * x + y + w - z, data[:a_c], data[:y_c], data[:b], data[:a_fc]));
data[:w_o] = (map((x,y) -> p.R * x + y, data[:a_p], data[:y_p]));
data[:w_y] = (map((x,y,z) -> p.R * x + y + z, data[:a_c], data[:y_c], data[:b]));

# return the asset (or wealth) of the stationary distribution
function getAsset(Fam::Array{Household})
  asset_o = [Fam[i].a_p for i in 1:length(Fam)]
  asset_y = Float64[]
  for h in Fam
    if h.t == 4
      push!(asset_y, h.a_c)
    end
  end
  asset   = Float64[]
  append!(asset,asset_y)
  append!(asset,asset_o)
  return (asset, asset_o, asset_y)
end

function getWealth(Fam::Array{Household})
  asset_o = [Fam[i].a_p + p.Y[Fam[i].yi_p,Fam[i].t] for i in 1:length(Fam)]
  asset_y = Float64[]
  for h in Fam
    push!(asset_y, h.a_c + p.Y[h.yi_c,h.t-2])
  end
  asset   = Float64[]
  append!(asset,asset_y)
  append!(asset,asset_o)
  return (asset, asset_o, asset_y)
end

# compute what the top percentage holds as a fraction of the total wealth
# e.g. for top 1%, do comp_top2(...,99)
function comp_top2(asset::Array{Float64}, perc::Int)
  total = sum(asset)
  q     = percentile(asset,perc)
  own   = 0.
  for a in asset
    if a > q
      own += a
    end
  end
  return own / total
end

# wealth histograms used in the thesis
asset_al, asset_o_al, asset_y_al = getWealth(fam_al)
asset_noal, asset_o_noal, asset_y_noal = getWealth(fam_noal);

plt1 = Plots.histogram(asset_noal, nbins = 5, title = "Aggregate", xlab = "Wealth",
                      label = "No altruism", grid = false, alpha = 0.3, color = :red,
                      yformatter = :scientific, legendfont = font(11), titlefont = font(14))
Plots.histogram!(asset_al, nbins = 20, label = "Altruism", alpha = 0.7, color = :darkred)

plt2 = Plots.histogram(asset_y_noal, nbins = 5, title = "Youngs", xlab = "Wealth",
                      label = "No altruism", grid = false, alpha = 0.3, color = :red,
                      yformatter = :scientific, legendfont = font(11), titlefont = font(14))
Plots.histogram!(asset_y_al, nbins = 20, label = "Altruism", alpha = 0.7, color = :darkred)

plt3 = Plots.histogram(asset_o_noal, nbins = 5, title = "Olds", xlab = "Wealth",
                      label = "No altruism", grid = false, alpha = 0.3, color = :red,
                      yformatter = :scientific, legendfont = font(11), titlefont = font(14))
Plots.histogram!(asset_o_al, nbins = 20, label = "Altruism", alpha = 0.7, color = :darkred)
l = @layout([a; b c ])
Plots.plot(plt1,plt2,plt3,layout = l)
