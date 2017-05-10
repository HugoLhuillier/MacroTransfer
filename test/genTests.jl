module genTests

include("./../src/GenModel/PFI.jl")
using Base.Test

p   = GenPFI.Param(size=10)
pol = GenPFI.Policies(p)
U   = GenPFI.Utility(p)

info("Starting testing on general model")

# savings and transfer guess used in the initialization phase
sav(a,b,y) = (y + b + p.R * a) / (1 + p.R)
tra(a,y)   = p.α * sav(a,0.,y)

@testset "re-usable functions" begin
  @testset "u and u'" begin
    @test U.du(3.) == 1/3
    p   = GenPFI.Param(size=10, ut = "CRRA", ɣ = 3.)
    U   = GenPFI.Utility(p)
    @test U.du(3.) == 3.^(-p.ɣ)
  end

  # this test will be included when the problem with AxisAlgorithm.jl will be resolved. see README.md
  @testset "interpolant" begin
    # can interplations with cubic splines match linear functions?
    # can it get that the our guesses are independent of some variables, and increasing in others?
    # warn("Turn these off before commiting for Travis")
    # A1, A2, B4 = GenPFI.interp(p,pol.A1,"A1"), GenPFI.interp(p,pol.A2,"A2"), GenPFI.interp(p,pol.B4,"B4")
    # @test A1[1,1][1.,3.] == A1[1,1][2.,3.]
    # @test A1[1,1][1.,2.] < A1[1,1][1.,3.]
    # @test pol.A1[10,10,1,1] - 1e-9 <= A1[1,1][p.a_grid[10], p.b_grid[10]] <= pol.A1[10,10,1,1] + 1e-9
    # @test B4[1,1][1.,3.] == B4[1,1][1.,2.]
    # @test B4[1,1][2.,3.] < B4[1,1][3.,3.]
    # @test pol.B4[10,10,1,1] - 1e-9 <= B4[1,1][p.a_grid[10], p.a_grid[10]] <= pol.B4[10,10,1,1] + 1e-9
    #
    # # test whether the extrapolation is exploding or not
    # @test_approx_eq A1[1,1][1,p.b_grid[end]+5] sav(0.,p.b_grid[end]+5,p.Y[1,1])
    # @test_approx_eq A2[1][p.a_grid[end]+5,1,0.] sav(p.a_grid[end]+5,0.,p.Y[1,2])
    # @test_approx_eq A2[1][0.,1,p.b_grid[end]+5] sav(0.,p.b_grid[end]+5,p.Y[1,2])
  end
end

@testset "initialization" begin
  # initial guess: savings (and transfers, when applicable) are increasing in own wealth,
  # transfer (if applicable) and income. they are indepedendent of the other agent's wealth,
  # and income. check for these
  # test also whether the dimensions are correct
  As, Bs, Ys = length(p.a_grid), length(p.b_grid), size(p.Y)[1]

  @testset "Age 1 policies" begin
    @test size(pol.A1)    == (As,Bs,Ys,Ys)
    @test pol.A1[1,:,:,:] == pol.A1[2,:,:,:];
    @test pol.A1[:,:,:,1] == pol.A1[:,:,:,1];
    @test all(x -> x == true, pol.A1[:,2,:,:] .> pol.A1[:,1,:,:]);
    @test all(x -> x == true, pol.A1[:,:,2,:] .> pol.A1[:,:,1,:]);
    @test pol.A1[1,end,1,1] == sav(0.,p.b_grid[end],p.Y[1,1])
  end

  @testset "Age 2 policies" begin
    @test size(pol.A2)      == (As,As,Bs,Ys)
    @test pol.A2[:,2,:,:]   == pol.A2[:,1,:,:];
    @test all(x -> x == true, pol.A2[:,:,2,:] .> pol.A2[:,:,1,:]);
    @test all(x -> x == true, pol.A2[:,:,:,2] .> pol.A2[:,:,:,1]);
    @test all(x -> x == true, pol.A2[2,:,:,:] .> pol.A2[1,:,:,:]);
    @test pol.A2[end,1,end,1,1] == sav(p.a_grid[end],p.b_grid[end],p.Y[1,2])
  end

  @testset "Age 3 policies" begin
    @test size(pol.A3)  == (As,Ys,Ys)
    @test size(pol.B3)  == (As,Ys,Ys)
    @test pol.A3[:,:,1] == pol.A3[:,:,2];
    @test all(x -> x == true, pol.A3[2,:,:] .> pol.A3[1,:,:]);
    @test all(x -> x == true, pol.A3[:,2,:] .> pol.A3[:,1,:]);
    @test pol.A3[end,1,1] == sav(p.a_grid[end],0.,p.Y[1,3])
    # approx because we previously multiply the matrix by p.α instead of the
    # the actual numbers, which lead to approximation at the 1e-9 level
    @test_approx_eq pol.B3[end,1,1] tra(p.a_grid[end],p.Y[1,3])
  end

  @testset "Age 4 policies" begin
    @test size(pol.A4)    == (As,As,Ys,Ys)
    @test size(pol.B4)    == (As,As,Ys,Ys)
    @test pol.A4[:,:,:,1] == pol.A4[:,:,:,2];
    @test pol.A4[:,2,:,:] == pol.A4[:,1,:,:];
    @test all(x -> x == true, pol.A4[2,:,:,:] .> pol.A4[1,:,:,:]);
    @test all(x -> x == true, pol.A4[:,:,2,:] .> pol.A4[:,:,1,:]);
    @test pol.A4[end,1,1,1] == sav(p.a_grid[end],0.,p.Y[1,4])
    @test_approx_eq pol.B4[end,1,1,1] tra(p.a_grid[end],p.Y[1,4])
  end

  @testset "Income profile" begin
    # want to check, for each age, top income > bottom income, and variance increases over time
    for i in 1:4
      @test p.Y[1,i] < p.Y[end,i]
      if i < 4
        @test var(p.Y[:,i]) < var(p.Y[:,i+1])
      end
    end
  end
end
# this test will be included when the problem with AxisAlgorithm.jl will be resolved. see README.md
# @testset "update" begin
#   @testset "Age 4 update" begin
#     A4_0, B4_0 = copy(pol.A4), copy(pol.B4)
#     GenPFI.old4!(p,U,pol)
#     # check that the mean has changed
#     @test mean(A4_0) != mean(pol.A4)
#     @test mean(B4_0) != mean(pol.B4)
#     # check that there is no savings and transfer below 0
#     @test all(x -> x == true, pol.A4 .>=0)
#     @test all(x -> x == true, pol.B4 .>=0)
#   end
#
#   @testset "Age 1 update" begin
#     A1_0 = copy(pol.A1)
#     GenPFI.young1!(p,U,pol)
#     # check that the mean has changed
#     @test mean(A1_0) != mean(pol.A1)
#     # check that there is no savings below 0
#     @test all(x -> x == true, pol.A1 .>=0)
#   end
#
#   @testset "Age 3 update" begin
#     A3_0, B3_0 = copy(pol.A3), copy(pol.B3)
#     GenPFI.old3!(p,U,pol)
#     # check that the mean has changed
#     @test mean(A3_0) != mean(pol.A3)
#     @test mean(B3_0) != mean(pol.B3)
#     # check that there is no savings and transfer below 0
#     @test all(x -> x == true, pol.A3 .>=0)
#     @test all(x -> x == true, pol.B3 .>=0)
#   end
# end

end
