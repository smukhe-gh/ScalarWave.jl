
#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
# Test summation-by-parts property for our derivative operator
# See Summation by Parts, Projections, and Stability. I Olsson 1995
#--------------------------------------------------------------------

D  = [-1.0  1.0  0.0  0.0;
      -0.5  0.0  0.5  0.0;
       0.0 -0.5  0.0  0.5;
       0.0  0.0 -1.0  1.0]

W  = [ 0.5  0.0  0.0  0.0;
       0.0  1.0  0.0  0.0;
       0.0  0.0  1.0  0.0;
       0.0  0.0  0.0  0.5]

B  = [-1.0  0.0  0.0  0.0;
       0.0  0.0  0.0  0.0;
       0.0  0.0  0.0  0.0;
       0.0  0.0  0.0  1.0] 

u = [x^2 for x in linspace(-1, 1, 4)]
v = [x^3 for x in linspace(-1, 1, 4)]
@test W*D + (W*D)' == B
@test u'*W*D*v + u'*(W*D)'*v == u'*B*v

# Now test the spectral discretization operators 
@testset "Summation-by-Parts" begin
    @testset "N = $N" for N in 1:6
        D = Float64[chebd(i,j, N) for i in 1:N+1, j in 1:N+1]
        W = diagm(Float64[chebw(i, N) for i in 1:N+1])
        B = zeros(N+1, N+1)
        B[1,1] = 1
        B[end, end] = -1
        u = [x^2 for x in chebgrid(N)]
        v = [x^3 for x in chebgrid(N)]
        @test u'*W*D*v + u'*(W*D)'*v ≈ u'*B*v atol=1e-15 
        @test_broken (W*D) + (W*D)' ≈ B atol=1e-15 
    end
end;
