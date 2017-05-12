# interpolate a policy function along the continuous states
# `pol::String` : the name policy to be interpolated
# `X::Array` : the policy function to be interpolated
function interp(p::Param, X::Array{Float64}, pol::String)
  # NOTE: cannot use linear interpolation if interpolant include NoInterp(), see issue #156
  # => hardcore the NoInterp() part of the interpolant
  Ys = size(p.Y)[1]
  if pol == "A1"
    F = Array{Any}(Ys,Ys)
    for i in 1:Ys
      for j in 1:Ys
        # itp    = interpolate(X[:,:,i,j], (BSpline(Cubic(Line())),
        #                       BSpline(Cubic(Line()))), OnGrid())
        itp    = interpolate(X[:,:,i,j], (BSpline(Quadratic(Line())),
                              BSpline(Quadratic(Line()))), OnGrid())
        # itp    = interpolate(X[:,:,i,j], (BSpline(Linear()),
        #                       BSpline(Linear())), OnGrid())
        itp    = scale(itp, p.a_grid, p.b_grid)
        F[i,j] = extrapolate(itp, Linear())
      end
    end
  elseif pol == "A2"
    F = Array{Any}(Ys)
    for i in 1:Ys
      # itp  = interpolate(X[:,:,i], (BSpline(Cubic(Line())),
      #                                 BSpline(Cubic(Line()))), OnGrid())
      itp  = interpolate(X[:,:,i], (BSpline(Quadratic(Line())),
                                      BSpline(Quadratic(Line()))), OnGrid())
      # itp  = interpolate(X[:,:,i], (BSpline(Linear()),
      #                                 BSpline(Linear())), OnGrid())
      itp  = scale(itp, p.a_grid, p.b_grid)
      F[i] = extrapolate(itp, Linear())
    end
  elseif (pol == "A3") || (pol == "B3")
    F = Array{Any}(Ys,Ys)
    for i in 1:Ys
      for j in 1:Ys
        # itp    = interpolate(X[:,i,j], (BSpline(Cubic(Line()))), OnGrid())
        itp    = interpolate(X[:,i,j], (BSpline(Quadratic(Line()))), OnGrid())
        # itp    = interpolate(X[:,i,j], (BSpline(Linear())), OnGrid())
        itp    = scale(itp, p.a_grid)
        F[i,j] = extrapolate(itp, Linear())
      end
    end
  elseif pol == "B4"
    F = Array{Any}(Ys,Ys)
    for i in 1:Ys
      for j in 1:Ys
        # itp    = interpolate(X[:,:,i,j], (BSpline(Cubic(Line())),
        #                                BSpline(Cubic(Line()))), OnGrid())
        itp    = interpolate(X[:,:,i,j], (BSpline(Quadratic(Line())),
                                       BSpline(Quadratic(Line()))), OnGrid())
        # itp    = interpolate(X[:,:,i,j], (BSpline(Linear()),
        #                                BSpline(Linear())), OnGrid())
        itp    = scale(itp, p.a_grid, p.a_grid)
        F[i,j] = extrapolate(itp, Linear())
      end
    end
  else
    error("Policies can only be A1, A2, A3, B3, A4, B4")
  end

  return F
end
