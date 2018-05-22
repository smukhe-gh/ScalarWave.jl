#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 05-2018
#--------------------------------------------------------------------

epsilon = 3e-2
for U in (pi/2 - epsilon)*chebgrid(20), V in (pi/2 - epsilon)*chebgrid(20)
    (r,t) = find_TR_from_UV(U, V)
    @test find_UV_from_TR(r,t)[1] ≈ U
    @test find_UV_from_TR(r,t)[2] ≈ V
end

params = Params((4,4), (-pi/2, pi/2), (-pi/2, pi/2), 1)
@test createmesh(params)[2,4][1] == (pi/2)*chebx(2,3)

