#==========================#
# VFI & PFI
#==========================#

# interpolate a policy function along one dimension (asset or transfer), and include
# a lookup table for the wage.
# `pol::Char` : `a` if want to interpolate the asset policy function, `b` for the transfer
# `x::AbstractMatrix` : Matrix of the policy function to be interpolated
function itp_pol(p::Param, pol::Char, x::AbstractMatrix)

  siz      = size(x)
  itp      = interpolate(x, (BSpline(Linear()), NoInterp()), OnGrid())
  if pol == 'a'
    # Interpolate the asset policy function, therefore along the transfer grid
    itp_s  = scale(itp, p.b_grid, 1:siz[2])
  elseif pol == 'b'
    # Interpolate the transfer policy function, therefore along the asset grid
    itp_s  = scale(itp, p.a_grid, 1:siz[2])
  end

  return itp_s
end

#==========================#
# EGM
#==========================#

# for a sorted vector, find the indeces of its first elements that brackets x
function findindex_bound(vect::Vector, x::Real)
  if x > maximum(vect)
    # NOTE: if x exceeds the max. of the vector, then returns the two highest indeces
    x0, x1 = Void, length(vect)
  else
    for i in 1:(length(vect)-1)
      if (vect[i] <= x) && (vect[i+1] > x)
        x0 = i
        x1 = i + 1
        break
      end
    end
  end
  return (x0, x1)
end

# linearly interpolate a function F, defined on X, evaluated at xi
function linear_itp(xi::Float64, X::Vector, F::Vector)
  x0, x1 = findindex_bound(X, xi)
  if x0 == Void
    x0     = x1 - 1
    f1i    = F[x1] + (xi - X[x1]) * (F[x1] - F[x0]) / (X[x1] - X[x0])
  else
    f1i    = F[x0] + (xi - X[x0]) * (F[x1] - F[x0]) / (X[x1] - X[x0])
  end
  return f1i
end

# interpolate the young's policy function in the transfer dimension, with a
# lookup table for the income, and return the derivative
function der_pol(p::Param, x::AbstractMatrix)
  siz      = size(x)
  itp      = interpolate(x, (BSpline(Linear()), NoInterp()), OnGrid())
  itp      = scale(itp, p.b_grid, 1:siz[2])
  d        = (xi, yi) -> gradient(itp, xi, yi)[1]
  return d
end
