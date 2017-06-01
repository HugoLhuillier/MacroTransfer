A2     = GenPFI.interp(p, pol.A2, "A2");
A3, B3 = GenPFI.interp(p, pol.A3, "A3"), GenPFI.interp(p, pol.B3, "B3");
A4, B4 = GenPFI.interp(p, pol.A4, "A4"), GenPFI.interp(p, pol.B4, "B4");

function getParam1(ai_fp,bi,yi_c,yi_p)
  return (p.a_grid[ai_fp], p.b_grid[bi], p.Y[yi_c,1], yi_c, yi_p)
end

function f1(x::Float64, a_fp::Float64, b::Float64, y_c::Float64, yi_c::Int, yi_p::Int)
  if x > y_c + b
    return 1e9
  else
    Eu2, Eu3 = 0., 0.
    for yi_fc in 1:length(p.Y[:,2])
      for yi_fp in 1:length(p.Y[:,4])
        # interpolate the parent 4 policies
        _b4 = B4[yi_fp,yi_fc][a_fp, x]
        _a4 = A4[yi_fp,yi_fc][a_fp, x]
        # ensure the non-negativity / corner of the Markov
        if _b4 < 1e-3
          _b4, _db = 0., 0.
        else
          _db  = gradient(B4[yi_fp,yi_fc], a_fp, x)[2]
        end
        if _a4 < 1e-3
          _a4, _da = 0., 0.
        else
          _da  = gradient(A4[yi_fp,yi_fc], a_fp, x)[2]
        end
        # own future policy
        _a2 = A2[yi_fc][x,_a4,_b4]
        _a2 = _a2 < 1e-3 ? 0. : _a2
        # expected period 2 utility, taking into account the probability of the incomes
        Eu2 += p.П[yi_c,yi_fc,2] * p.П[yi_p,yi_fp,4] * U.du(p.R * x + p.Y[yi_fc,2] + _b4 - _a2) * (p.R + _db)

        # compute period 3 expected utility
        for yi_ffp in 1:length(p.Y[:,3])
          for yi_ffc in 1:length(p.Y[:,1])
            # own choices
            _a3 = A3[yi_ffp,yi_ffc][_a2 + _a4]
            _b3 = B3[yi_ffp,yi_ffc][_a2 + _a4]
            _a3 = _a3 < 1e-3 ? 0. : _a3
            _b3 = _b3 < 1e-3 ? 0. : _b3

            # expected period 2 utility, taking into account the probability of the incomes
            Eu3 += U.du(p.R * (_a2 + _a4) + p.Y[yi_ffp,3] - _b3 - _a3) * _da * p.R * p.β *
                  p.П[yi_c,yi_fc,2] * p.П[yi_p,yi_fp,4] * p.П[yi_fc,yi_ffp,3] * p.П[yi_fc,yi_ffc,1]
          end
        end

      end
    end
    return U.du(y_c + b - x) - p.β * (Eu2 + Eu3)
  end
end

function foc1(par)
  w   = par[2] + par[3]
  X   = linspace(0.,w - w/5,100)
  fig = figure("FOC 1",figsize=(10,10))
  ax  = fig[:add_subplot](1,1,1)
  ax[:plot](X,f1.(X,par...), color = "red")
  ax[:plot](X, zeros(100), color = "black", linestyle="--")
  xlabel(L"$\tilde{a}'$")
  title("FOC #1")
end
