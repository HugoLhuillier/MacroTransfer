module ToyPlot

include("./../ToyModel/ToyModel.jl")
using PyPlot

function cf(;save::Bool=false, isPFI::Bool=true)

  p              = ToyModel.Param(grid_size = 20, max_a=8)
  y_p, y_r       = p.inc
  a_grid, b_grid = p.a_grid, p.b_grid
  beta, al, R    = p.beta, p.alpha, p.R
  if isPFI
    sol          = ToyModel.PFI(20,p)
  else
    sol          = ToyModel.EGM(20,p)
  end

  # build the Markov strategy
  a_p(b::Float64) = beta * (y_p + b) / (1 + beta)
  b_p(a::Float64) = 0.
  function b_r(a::Float64)
      if al * R * a > y_p
          return (al * R * a - y_p) / (1 + al)
      else
          return 0
      end
  end
  function a_r(b::Float64)
      w     = b + y_r
      if w <= (1 + beta) * y_p / (al * beta * R)
          return beta * w / (1 + beta)
      else
          return (beta * R * w * (1 + al) - y_p) / (R * (1 + beta * (1 + al)))
      end
  end

  # Equilibrium level
  as_r = [(y_r * beta * R * (1 + al) - y_p) / (R * (1 + beta * (1 + al)))]
  as_p = [al * beta * beta * (y_p + R * y_r) / ((1 + beta) * 1 + beta * (1 + al))]
  bs_r = [0.]
  bs_p = [(al * beta * R * y_r - (1 + beta) * y_p) / (1 + beta * (1 + al))]

  fig = figure("cf_vs_num",figsize=(10,10))
  subplots_adjust(wspace=0.4, hspace=0.4)
  ax  = fig[:add_subplot](2,2,1)
  ax[:scatter](bs_r, a_r(bs_r[1]), label = "Eq.", color = "darkred", zorder = 3)
  ax[:plot](b_grid, sol[1][:,2], label = "Num.", color = "indianred", zorder = 1)
  ax[:plot](b_grid, a_r.(b_grid), linestyle=":", label = "Analytical", zorder = 2, color = "mistyrose")
  xlabel(L"$b$")
  title("Savings rich")
  legend()
  ax  = fig[:add_subplot](2,2,2)
  ax[:plot](b_grid, sol[1][:,1], label = "Num.", color = "indianred", zorder = 1)
  ax[:plot](b_grid, a_p.(b_grid), linestyle=":", label = "Analytical", color = "mistyrose", zorder = 2)
  ax[:scatter](bs_p, a_p(bs_p[1]), label = "Eq.", color = "darkred", zorder = 3)
  xlabel(L"$b$")
  title("Savings poor")
  legend()
  ax  = fig[:add_subplot](2,2,3)
  ax[:plot](a_grid, sol[2][:,2], label = "Num.", color = "indianred", zorder = 1)
  ax[:plot](a_grid, b_p.(b_grid), linestyle=":", label = "Analytical", color = "mistyrose", zorder = 2)
  ax[:scatter](as_p, b_p(as_p[1]), label = "Eq.", color = "darkred", zorder = 3)
  xlabel(L"$a$")
  title("Transfer to rich")
  yticks([0.])
  legend()
  ax  = fig[:add_subplot](2,2,4)
  ax[:plot](a_grid, sol[2][:,1], label = "Num.", color = "indianred", zorder = 1)
  ax[:plot](a_grid, b_r.(a_grid), linestyle=":", label = "Analytical", color = "mistyrose", zorder = 2)
  ax[:scatter](as_r, b_r(as_r[1]), label = "Eq.", color = "darkred", zorder = 3)
  xlabel(L"$a$")
  title("Transfer to poor")
  legend()

  if save
    PyPlot.savefig("num_vs_cf.pdf")
  end
end

function accuracyVFI()
  p     = ToyModel.Param(grid_size=100, max_a = 8.)
  sol_v = ToyModel.VFI(100,p)
  sol_p = ToyModel.PFI(100,p)

  fig = figure("cf_vs_num",figsize=(10,10))
  ax  = fig[:add_subplot](1,2,1)
  ax[:plot](p.b_grid, sol_v[1][:,1], linestyle=":", color = "red")
  ax[:plot](p.b_grid, sol_v[1][:,2], linestyle=":", color = "darkred")
  ax[:plot](p.b_grid, sol_p[1][:,1], color = "red")
  ax[:plot](p.b_grid, sol_p[1][:,2], color = "darkred")
  xlabel(L"$b$")
  legend(["VFI, poor", "VFI, rich", "Ana., poor", "Ana. rich"])
  title("Savings")
  legend()

  ax  = fig[:add_subplot](1,2,2)
  ax[:plot](p.a_grid, sol_v[2][:,1], linestyle=":", color = "red")
  ax[:plot](p.a_grid, sol_v[2][:,2], linestyle=":", color = "darkred")
  ax[:plot](p.a_grid, sol_p[2][:,1], color = "red")
  ax[:plot](p.a_grid, sol_p[2][:,2], color = "darkred")
  xlabel(L"$a$")
  legend(["VFI, poor", "VFI, rich", "Ana., poor", "Ana. rich"])
  title("Transfer")
  legend()
end

function effiency(;grid_max=500)
  max_g = [20:100:grid_max;grid_max]
  vfi_t = zeros(length(max_g))
  pfi_t = zeros(length(max_g))
  egm_t = zeros(length(max_g))

  # pre-compile the algorithm
  p = ToyModel.Param(grid_size=10)
  ToyModel.PFI(30,p)
  ToyModel.EGM(30,p)
  ToyModel.VFI(100,p)

  for (i, a) in enumerate(max_g)
    println("Currently @ $i out of $(length(max_g))")
    p = ToyModel.Param(grid_size=a)
    pfi_t[i] = @elapsed ToyModel.PFI(30,p);
    egm_t[i] = @elapsed ToyModel.EGM(30,p);

    t = @elapsed sol = ToyModel.VFI(200,p);
    if typeof(sol) <: Tuple
      vfi_t[i] = t
    else
      vfi_t[i] = NaN
    end
  end

  fig = figure("effiency",figsize=(10,10))
  ax  = fig[:add_subplot](1,1,1)
  ax[:scatter](max_g,vfi_t,label="VFI")
  ax[:scatter](max_g,pfi_t,label="PFI")
  ax[:scatter](max_g,egm_t,label="EGM")
  xlabel("Grid size")
  ylabel("Time (s)")
  legend()
end

end
