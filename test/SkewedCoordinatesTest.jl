#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 03-2019
# Test wave equation on Minkowski in 3D
#--------------------------------------------------------------------

using LinearAlgebra
struct W end

SW = GaussLobatto{W, 15, 1.0, -1.0} 
SV = GaussLobatto{V, 15, 1.0, -1.0}
SU = GaussLobatto{U, 15, 1.0, -1.0}
SUVW  = ProductSpace{SW, SV, SW}
DW, DV, DU = derivative(SUVW)

D  = [DU, DV, DW]
B  = boundary(Null, SUVW)
u0 = Field(SUVW, (U,V,W) -> exp(-((U - ((maximum(SU) + minimum(SU))/2))^2)/0.1)
                          + exp(-((V - ((maximum(SV) + minimum(SV))/2))^2)/0.1)
                          + exp(-((W - ((maximum(SW) + minimum(SW))/2))^2)/0.1))

u0 = Field(SUVW, (U,V,W) -> exp(-U^2/0.1) + exp(-W^2/0.1)  + exp(-V^2/0.1)) 
τ  = Field(SUVW, (U,V,W) -> (π/4)*cospi(U/2)*cospi(V/2)*cospi(W/2)) 

# g = Metric{_dd, 3}([0, -1, -1,
                        # 0, -1,
                            # 0])

g = Metric{_uu,3}([-0*τ, -cos(τ) - sin(τ),  -cos(τ) + sin(τ),
                         -2*cos(τ)*sin(τ),         -cos(2*τ),
                                                    sin(2*τ)])
                  

# g = Metric{_dd, 3}([1,  0,    0,
                     # -1/3,  2/3,
                            # 1/6])

L = sum(g[i,j]*D[i]*D[j] for i in 1:3, j in 1:3)
using LinearAlgebra
@show cond(reshape(L ⊙ B))
u = solve(L ⊙ B, B*u0)

@show maximum(abs.(u.value[1, :, :] - u.value[end, :, :]))
@show maximum(abs.(u.value[:, 1, :] - u.value[:, end, :]))
@show maximum(abs.(u.value[:, :, 1] - u.value[:, :, end]))

# Before moving ahead, check if the boundary is preserved. 
# @test_broken B*u ≈ B*u0 
@show cond(reshape(L⊙B))
# @show cond(reshape(L⊙B + B*DU + B*DV + B*DW))
@show maximum(abs(B*u - B*u0))
writevtk(u, "../output/wave-solution")
writevtk(B*u, "../output/wave-solution-boundary")

