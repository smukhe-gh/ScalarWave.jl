#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
# Basis transformation Test
#--------------------------------------------------------------------

@testset "coordinate transformation" begin
    for i in 1:3, j in 1:3
        u = -abs(4*randn())
        v =  abs(4*randn())

        t = find_t_of_UV(u,v,1.0)
        r = find_r_of_UV(u,v,1.0)
        u_of_tr = find_U_of_tr(t,r,1.0) 
        v_of_tr = find_V_of_tr(t,r,1.0) 

        @test u ≈ u_of_tr
        @test v ≈ v_of_tr
    end
end
