#============================#
#       STAGE GAME
#============================#

# stores the information relevant to the family
type Household
  a_p::Float64
  a_c::Float64
  b::Float64
  yi_p::Int64
  yi_c::Int64
  t::Int64

  # create a new instance
  function Household(t::Int64, p::GenPFI.Param)
    yi_p = rand(1:2)
    yi_c = rand(1:2)
    a_p  = p.a_grid[rand(1:length(p.a_grid))]
    return new(a_p,0.,0.,yi_p,yi_c,t)
  end
end

# simulate the income shock
function sim_shock!(p::GenPFI.Param, h::Household)
  P_o = DiscreteRV(vec(p.П[h.yi_p,:,h.t]))
  P_y = DiscreteRV(vec(p.П[h.yi_c,:,h.t-2]))
  h.yi_p = QuantEcon.draw(P_o)
  h.yi_c = QuantEcon.draw(P_y)
  return nothing
end

# given the income, return the best response
function stage_game!(p::GenPFI.Param, h::Household, itp::Dict)
  if h.t == 3
    A3, B3, A1 = itp["A3"], itp["B3"], itp["A1"]
    _a3   = A3[h.yi_p,h.yi_c][h.a_p]
    h.a_p = p.R * (_a3 < 0. ? 0. : _a3)
    _b3   = B3[h.yi_p,h.yi_c][h.a_p]
    h.b   = _b3 < 0. ? 0. : _b3
    _a1   = A1[h.yi_c,h.yi_p][h.a_p, h.b]
    h.a_c = p.R * (_a1 < 0. ? 0. : _a1)
  elseif h.t == 4
    B4, A2 = itp["B4"], itp["A2"]
    _b4   = B4[h.yi_p,h.yi_c][h.a_p,h.a_c]
    h.b   = _b4 < 0. ? 0. : _b4
    _a2   = A2[h.yi_c][h.a_c,h.b]
    h.a_p = p.R * (_a2 < 0. ? 0. : _a2)
    h.a_c = 0.
  else
    error("t can only be 3 or 4")
  end
  return nothing
end

# update the age profile
function age!(p::GenPFI.Param, h::Household)
  if h.t == 3
    h.t = 4
  elseif h.t == 4
    h.t = 3
  else
    error("t can only be 3 or 4")
  end
  return nothing
end

# simulate the stage game, for a given household
function update_h!(p::GenPFI.Param, h::Household, itp::Dict)
  sim_shock!(p,h)
  stage_game!(p,h,itp)
  age!(p,h)
end

#============================#
#       MEASURE
#============================#

# compute the moments of the measure λ. use as moments
# mean, standard deviation, skewness and kurtosis, as well as the 20, 40, 60 and 80 quantile
# mean(x), std(x), skewness(x), kurtosis(x), nquantile(x,5)[2:5]
function moments(λ::Array{Household})
  pop_o = [λ[i].a_p for i in 1:length(λ)]
  pop_y = [λ[i].a_c for i in 1:length(λ)]
  m1    = mean(pop_o)
  m2    = std(pop_o)
  m3    = skewness(pop_o)
  m4    = kurtosis(pop_o)
  m5, m6, m7, m8 = nquantile(pop_o,5)[2:5]

  m9   = mean(pop_y)
  m10  = std(pop_y)
  m11  = skewness(pop_y)
  m12  = kurtosis(pop_y)
  m13, m14, m15, m16 = nquantile(pop_y,5)[2:5]

  return Dict("mean_o" => m1, "std_o" => m2, "skew_o" => m3, "kurt_o" => m4,
              "q20_o" => m5, "q40_o" => m6, "q60_o" => m7, "q80_o" => m8,
              "mean_y" => m9, "std_y" => m10, "skew_y" => m11, "kurt_y" => m12,
              "q20_y" => m13, "q40_y" => m14, "q60_y" => m15, "q80_y" => m16)
end

# store information relevant to the measure: the individuals, and the moments of the distribution
type Measure
  λ::Array{Household}
  mom::Dict

  function Measure(Families::Array{Household}, p::GenPFI.Param, t::Int)
    λ = Household[]
    for h in Families
      if h.t == t
        push!(λ, h)
      end
    end
    mom = moments(λ)
    return new(λ, mom)
  end
end

function update_meas!(Families::Array{Household},m3::Measure,m4::Measure)
  λ3, λ4  = Household[], Household[]
  for h in Families
    if h.t == 3
      push!(λ3,h)
    else
      push!(λ4,h)
    end
  end
  m3.λ = λ3
  m4.λ = λ4
  return nothing
end

# check the convergence of the distribution, based on some moments
function check_conv(m3::Measure, m4::Measure, tol::Float64)
  d3 = moments(m3.λ)
  d4 = moments(m4.λ)

  err = 0
  for k in keys(d3)
    if abs(d3[k] - m3.mom[k]) > err
      err = abs(d3[k] - m3.mom[k])
    end
  end
  for k in keys(d4)
    if abs(d4[k] - m4.mom[k]) > err
      err = abs(d4[k] - m4.mom[k])
    end
  end

  m3.mom = d3
  m4.mom = d4
  return err
end

#============================#
#       SIMULATION
#============================#

# build the initial Array of households
function init(numbFam::Int)
  # init the families
  Families = Array{Household}(numbFam)
  for i in 1:numbFam
    if i <= Int(numbFam/2)
      Families[i] = Household(3,p)
    else
      Families[i] = Household(4,p)
    end
  end
  return Families
end

# main function
function stationarity(p::GenPFI.Param, pol::GenPFI.Policies, numbFam::Int;
                      tol::Float64=1e-3, maxIter::Int=100)

  A1, A2, A3 = GenPFI.interp_lin(p, pol.A1, "A1"), GenPFI.interp_lin(p, pol.A2, "A2"), GenPFI.interp_lin(p, pol.A3, "A3");
  B3, B4     = GenPFI.interp_lin(p, pol.B3, "B3"), GenPFI.interp_lin(p, pol.B4, "B4");
  itp        = Dict("A1" => A1, "A2" => A2, "A3" => A3, "B3" => B3, "B4" => B4);

  # need to initialize the distribution
  Families = init(numbFam)
  m3, m4   = Measure(Families,p,3), Measure(Families,p,4)

  for iter in 1:maxIter
    for h in Families
      update_h!(p, h, itp)
    end
    update_meas!(Families,m3,m4)
    err      = check_conv(m3,m4,tol)
    info("@ iter $iter, error = $err")
    if err < tol
      return Families
    elseif iter == maxIter
      error("Did not converged in $maxIter iterations")
    end
  end
end

#============================#
#       PANEL SIMULATION
#============================#

function sim!(h::Household, i::Int, p::GenPFI.Param, sim::DataFrame, itp::Dict, t::Int)
  # store the state
  a_p = h.a_p
  a_c = h.a_c
  g   = h.t
  a_fp, a_fc, b = 0., 0., 0.
  # draw the income
  P_o = DiscreteRV(vec(p.П[h.yi_p,:,h.t]))
  P_y = DiscreteRV(vec(p.П[h.yi_c,:,h.t-2]))
  yi_p = QuantEcon.draw(P_o)
  yi_c = QuantEcon.draw(P_y)
  if h.t == 4
    var_p = p.Y[yi_p,4] - p.Y[h.yi_p,3]
    var_c = p.Y[yi_c,2] - p.Y[h.yi_c,1]
  else
    var_p = NaN
    var_c = NaN
  end
  h.yi_p = yi_p
  h.yi_c = yi_c
  y_p, y_c = p.Y[h.yi_p,g], p.Y[h.yi_c,g-2]
  # play the stage game
  if h.t == 3
    A3, B3, A1 = itp["A3"], itp["B3"], itp["A1"]
    a_fp  = A3[h.yi_p,h.yi_c][h.a_p]
    h.a_p = p.R * (a_fp < 0. ? 0. : a_fp)
    b     = B3[h.yi_p,h.yi_c][h.a_p]
    h.b   = b < 0. ? 0. : b
    a_fc  = A1[h.yi_c,h.yi_p][h.a_p, h.b]
    h.a_c = p.R * (a_fc < 0. ? 0. : a_fc)
  else h.t == 4
    B4, A2 = itp["B4"], itp["A2"]
    b     = B4[h.yi_p,h.yi_c][h.a_p,h.a_c]
    h.b   = b < 0. ? 0. : b
    a_fc  = A2[h.yi_c][h.a_c,h.b]
    h.a_p = p.R * (a_fc < 0. ? 0. : a_fc)
    a_fp  = 0.
    h.a_c = 0.
  end
  # then age for the simulation process
  if h.t == 3
    h.t = 4
  else h.t == 4
    h.t = 3
  end
  push!(sim, [a_p, a_c, a_fp, a_fc, b, y_p, var_p, y_c, var_c, g, i, t])
end

function simulate_pan(T::Int, Families::Array{Household}, p::GenPFI.Param, pol::GenPFI.Policies)
  sim = DataFrame(a_p = Float64[], a_c = Float64[], a_fp = Float64[], a_fc = Float64[],
                  b = Float64[], y_p = Float64[],
                  var_p = Float64[], y_c = Float64[],
                  var_c = Float64[], g = Int[], i = Int[], t = Int[])
  fam = sample(Families, 100_000)
  A1, A2, A3 = GenPFI.interp_lin(p, pol.A1, "A1"), GenPFI.interp_lin(p, pol.A2, "A2"), GenPFI.interp_lin(p, pol.A3, "A3");
  B3, B4     = GenPFI.interp_lin(p, pol.B3, "B3"), GenPFI.interp_lin(p, pol.B4, "B4");
  itp        = Dict("A1" => A1, "A2" => A2, "A3" => A3, "B3" => B3, "B4" => B4);
  for t in 1:T
    for (i,h) in enumerate(fam)
      sim!(h, i, p, sim, itp, t)
    end
  end
  return sim
end
