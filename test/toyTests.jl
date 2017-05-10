module toyTests

using Base.Test
include("./../src/ToyModel/ToyModel.jl")

# What we need to tests
# 1. There is indeed some update
# 2. The results correspond to the analical solution
# 3. The algorithm converges

info("Starting testing on toy model")
p    = ToyModel.Param(grid_size=20)

@testset "re-usable functions" begin
  @testset "u, u' and inverse u'" begin
    # same functions for PFI, VFI and EGM
    p   = ToyModel.Param(grid_size=20)
    u   = ToyModel.pick_u(p)
    du  = ToyModel.pick_du(p)
    idu = ToyModel.pick_idu(p)
    @test u(3.)   == log(3.)
    @test du(3.)  == 1/3
    @test idu(3.) == 1/3

    p   = ToyModel.Param(grid_size=20, ut = "CRRA", gamma = 2.)
    u   = ToyModel.pick_u(p)
    du  = ToyModel.pick_du(p)
    idu = ToyModel.pick_idu(p)
    @test u(3.)   == (3^(1-p.gamma) - 1) / (1 - p.gamma)
    @test du(3.)  == 3.^(-p.gamma)
    @test idu(3.) == (1/3)^(1/p.gamma)
  end

  @testset "interpolant" begin
    # only for PFI and VFI
    A1, B1 = ToyModel.init_P(p)

    A1_itp, B1_itp = ToyModel.itp_pol(p,'a',A1), ToyModel.itp_pol(p,'b',B1)
    @test_approx_eq A1[1,1] A1_itp[p.b_grid[1],1]
    for i in 1:5:20
      @test A1[i,1] - 1e-9 <= A1_itp[p.b_grid[i],1] <= A1[i,1] + 1e-9
      @test A1[i,2] - 1e-9 <= A1_itp[p.b_grid[i],2] <= A1[i,2] + 1e-9
      @test B1[i,1] - 1e-9 <= B1_itp[p.a_grid[i],1] <= B1[i,1] + 1e-9
      @test B1[i,2] - 1e-9 <= B1_itp[p.a_grid[i],2] <= B1[i,2] + 1e-9
    end

    X = collect(1:1:20)
    Y = [1;1;3;5;7]
    @test (3,4)            == ToyModel.findindex_bound(X,3)
    @test (Void,length(X)) == ToyModel.findindex_bound(X,22)
    @test (2,3)            == ToyModel.findindex_bound(Y,1)

    # linear function => linear interpolation should get it right
    f(x) = 2x + 1
    xx   = linspace(0.,20.,8)
    F    = f.(xx)
    # test both inside and outside the grid
    @test f(2.) == ToyModel.linear_itp(2.,collect(xx),F)
    @test f(30.) == ToyModel.linear_itp(30.,collect(xx),F)
  end
end

@testset "initialization" begin
  # PFI
  A, B = ToyModel.init_P(p)
  @test size(A) == (length(p.b_grid),length(p.inc))
  @test size(B) == (length(p.a_grid),length(p.inc))
  # VFI
  V0, A, B = ToyModel.init_V(p)
  @test size(A)  == (length(p.b_grid),length(p.inc))
  @test size(B)  == (length(p.a_grid),length(p.inc))
  @test size(V0) == (length(p.a_grid),length(p.inc))
  # EGM
  A, B = ToyModel.init_E(p)
  @test size(A) == (length(p.b_grid),length(p.inc))
  @test size(B) == (length(p.a_grid),length(p.inc))
end

@testset "update" begin
  @testset "PFI" begin
    A1, B1 = ToyModel.init_P(p)
    A2     = copy(A1)
    ToyModel.young_P!(B1,A2,p)
    B2 = ToyModel.old_P(B1,A2,p)
    # check that indeed, the algorithm updated the solutions
    # note that we cannot do pointwise checking, as some of the cells might not have moved
    @test mean(A1[:,1]) != mean(A2[:,1])
    @test mean(A1[:,2]) != mean(A2[:,2])
    @test mean(B1[:,1]) != mean(B2[:,1])
    @test mean(B1[:,2]) != mean(B2[:,2])
  end

  @testset "VFI" begin
    V0_1, A1, B1 = ToyModel.init_V(p)
    A2, B2       = copy(A1), copy(B1)
    ToyModel.young_V!(V0_1,A2,p)
    V0_2 = ToyModel.old_V!(V0_1,A2,B2,p)
    @test mean(A1[:,1])   != mean(A2[:,1])
    @test mean(A1[:,2])   != mean(A2[:,2])
    @test mean(B1[:,1])   != mean(B2[:,1])
    @test mean(B1[:,2])   != mean(B2[:,2])
    @test mean(V0_1[:,1]) != mean(V0_2[:,1])
    @test mean(V0_1[:,2]) != mean(V0_2[:,2])
  end

  @testset "EGM" begin
    CY0, CO0 = ToyModel.init_E(p)
    CY1 = copy(CY0)
    ToyModel.young_E!(p,CY1,CO0)
    CO1 = ToyModel.old_E(p,CY1,CO0)
    @test mean(CY0[:,1]) != mean(CY1[:,1])
    @test mean(CY0[:,2]) != mean(CY1[:,2])
    @test mean(CO0[:,1]) != mean(CO1[:,1])
    @test mean(CO0[:,2]) != mean(CO1[:,2])
  end
end

@testset "results" begin
  @testset "VFI" begin
    # test of convergence
    sol = ToyModel.VFI(100,p)
    # if did not converge, then return Void
    @test typeof(sol) <: Tuple{Array,Array}
  end

  p              = ToyModel.Param(grid_size = 20, max_a=8)
  y_p, y_r       = p.inc
  a_grid, b_grid = p.a_grid, p.b_grid
  beta, al, R    = p.beta, p.alpha, p.R

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

  # point at which to evaluate the solution
  a = [1;5;10;15;20]
  b = [1;5;10;15;20]

  @testset "PFI" begin
    # test of convergence
    sol = ToyModel.PFI(20,p)
    # if did not converge, then return Void
    @test typeof(sol) <: Tuple{Array,Array,Array,Array}
    # then test if solutions are equal to the analytical ones
    for i in 1:length(a)
      @test a_p(p.b_grid[i]) - 1e-6 <= sol[1][i,1] <= a_p(p.b_grid[i]) + 1e-6
      @test a_r(p.b_grid[i]) - 1e-6 <= sol[1][i,2] <= a_r(p.b_grid[i]) + 1e-6
      @test b_r(p.a_grid[i]) - 1e-6 <= sol[2][i,1] <= b_r(p.a_grid[i]) + 1e-6
      @test b_p(p.a_grid[i]) - 1e-6 <= sol[2][i,2] <= b_p(p.a_grid[i]) + 1e-6
    end
  end

  @testset "EGM" begin
    # convergence test
    sol = ToyModel.EGM(20,p)
    @test typeof(sol) <: Tuple{Array,Array,Array,Array}
    # then test if solutions are equal to the analytical ones
    for i in 1:length(a)
        @test a_p(p.b_grid[i]) - 1e-6 <= sol[1][i,1] <= a_p(p.b_grid[i]) + 1e-6
        @test a_r(p.b_grid[i]) - 1e-6 <= sol[1][i,2] <= a_r(p.b_grid[i]) + 1e-6
        @test b_r(p.a_grid[i]) - 1e-6 <= sol[2][i,1] <= b_r(p.a_grid[i]) + 1e-6
        @test b_p(p.a_grid[i]) - 1e-6 <= sol[2][i,2] <= b_p(p.a_grid[i]) + 1e-6
    end
  end
end

end
