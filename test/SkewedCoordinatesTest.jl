#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 03-2019
# Test wave equation on Minkowski in 3D
#--------------------------------------------------------------------

using LinearAlgebra
struct W end

SW = GaussLobatto{W, 18, 1.0, -1.0} 
SV = GaussLobatto{V, 18, 1.0, -1.0}
SU = GaussLobatto{U, 18, 1.0, -1.0}
SUVW  = ProductSpace{SW, SV, SW}
DW, DV, DU = derivative(SUVW)

D  = [DU, DV, DW]
B  = boundary(Null, SUVW)
u0 = Field(SUVW, (U,V,W) -> exp(-((U - ((maximum(SU) + minimum(SU))/2))^2)/0.1)
                          + exp(-((V - ((maximum(SV) + minimum(SV))/2))^2)/0.1)
                          + exp(-((W - ((maximum(SW) + minimum(SW))/2))^2)/0.1))

# Squished cube metric
# g = Metric{_dd, 3}([1/1, 0/1,   0/1,
                        # -1/3,   4/3,
                                # 2/3])

g = Metric{_dd, 3}([0, -1, -1,
                        0, -1,
                            0])

L = sum(g[i,j]*D[i]*D[j] for i in 1:3, j in 1:3)
u = solve(L ⊙ B, B*u0)

# Before moving ahead, check if the boundary is preserved. 
# @test_broken B*u ≈ B*u0 
@show cond(reshape(L⊙B))
# @show cond(reshape(L⊙B + B*DU + B*DV + B*DW))
@show maximum(abs(B*u - B*u0))
writevtk(u, "../output/wave-solution")
writevtk(B*u, "../output/wave-solution-boundary")

