"""
```
Policies(p::Param)
```

Constructor of the Policies object, storing the policy functions. In addition, it
initializes the policy functions.

#### State spaces

- Age 1: parent's asset x transfer x child's income x parent's income
- Age 2: child's asset x parent's asset x transfer x child's income
- Age 3: parent's asset x parent's income x child's income
- Age 4: parent's asset x child's asset x parent's income x child's income

*Initial guess*:

- Savings are linearly increasing in weatlh
- Transfers = savings adjusted by the altruistic parameter
"""
type Policies
  A1::Array{Float64}
  A2::Array{Float64}
  B3::Array{Float64}
  A3::Array{Float64}
  B4::Array{Float64}
  A4::Array{Float64}

  function Policies(p::Param)
    Ys = size(p.Y)[1]
    As = length(p.a_grid)
    Bs = length(p.b_grid)

    A1 = zeros(As,Bs,Ys,Ys)
    A2 = zeros(As,As,Bs,Ys)
    A3 = zeros(As,Ys,Ys)
    B3 = zeros(As,Ys,Ys)
    A4 = zeros(As,As,Ys,Ys)
    B4 = zeros(As,As,Ys,Ys)

    for (ai,a) in enumerate(p.a_grid)
      for (bi,b) in enumerate(p.b_grid)
        A1[:,bi,:,:]    = repeat(reshape(repeat((p.Y[:,1] + b) / (1 + p.R), inner=[As]),
                                        As, Ys, 1, 1),
                                outer = [1,1,1,Ys])
        A2[ai,:,bi,:] = reshape(repeat((p.Y[:,2] + b + p.R * a) / (1 + p.R), inner = [As]),
                                        As,Ys,1)
        A3[ai,:,:]      = reshape(repeat((p.Y[:,3] + p.R * a) / (1 + p.R), outer = [Ys]),
                                  Ys,Ys,1)
        B3[ai,:,:]      = p.α .* A3[ai,:,:]
        A4[ai,:,:,:]    = repeat(reshape(repeat((p.Y[:,4] + p.R * a) / (1 + p.R), inner = [As]),
                                          As,Ys,1),
                                outer = [1,1,Ys])
        B4[ai,:,:,:]    = p.α .* A4[ai,:,:,:]
      end
    end

    return new(A1,A2,B3,A3,B4,A4)
  end
end
