#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 03-2019
# Test wave equation on Minkowski in 3D
#--------------------------------------------------------------------

using LinearAlgebra
struct W end

SW = GaussLobatto{W, 8,  -0.75, -1.0} 
SV = GaussLobatto{V, 8,  -0.75, -1.0}
SU = GaussLobatto{U, 8,  -0.75, -1.0}
SUVW  = ProductSpace{SW, SV, SW}
DW, DV, DU = derivative(SUVW)

D = [DU, DV, DW]
B = boundary(Null, SUVW)
u = Field(SUVW, (U,V,W) -> exp(-((U - ((maximum(SU) + minimum(SU))/2))^2)/0.01)
                         + exp(-((V - ((maximum(SV) + minimum(SV))/2))^2)/0.01)
                         + exp(-((W - ((maximum(SW) + minimum(SW))/2))^2)/0.01))

b = B*u

writevtk(u, "bulk")
writevtk(b, "boundary")

guu = Field(SUVW, (U,V,W)->(1/3) - (1/sqrt(3)))
guv = Field(SUVW, (U,V,W)->(1/3)) 
guw = Field(SUVW, (U,V,W)->(-2/3)*sqrt(2 + sqrt(3)))
gvv = Field(SUVW, (U,V,W)->(1/3) + (1/sqrt(3)))
gvw = Field(SUVW, (U,V,W)->(sqrt(2)/3)*(-1 + sqrt(3)))
gww = Field(SUVW, (U,V,W)->(2/3)) 

g = Metric{_dd, 3}([guu, guv, guw,
                         gvv, gvw,
                              gww])

L = sum(g[i,j]*D[i]*D[j] for i in 1:3, j in 1:3)
s = prod(size(SUVW))^(2.0)^(3.0)
sinv = 1/prod(size(SUVW))^(2.0)^(3.0)
ev = eigvals(vec(L+B))

# @show maximum(abs.(eigvals(vec(L))))
# @show maximum(abs.(ev))
# @show minimum(abs.(ev))

@show cond(vec(L+B))
@show cond(vec(L + s*B))
@show cond(vec(L + sinv*B))
@show cond(vec(s*L + B))
@show cond(vec(sinv*L + B))

# using PyPlot
# plot(abs.(ev))
# show()

@show maximum(u.value)
u = solve(L + B, b)
@show maximum(u.value)
writevtk(u, "solution")
writevtk(B*u, "boundary-solution")
