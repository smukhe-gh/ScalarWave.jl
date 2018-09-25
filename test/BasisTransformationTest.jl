#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
# Basis transformation Test
#--------------------------------------------------------------------

struct M end
ϕ = Field(GaussLobatto{M, 9}, x->x^5 + 2)
ψ = basistransform(ϕ, Chebyshev{M, 9})
λ = basistransform(ψ, GaussLobatto{M, 9})

x = Field(GaussLobatto{M, 200}, x->x)
Φ = interpolate(ϕ, x)

#=
using PyPlot
plot(x.value, Φ.value)
show()
=#

struct N end
ω  = Field(ProductSpace{GaussLobatto{M,9}, GaussLobatto{N,12}}, (x,y)-> x^2 + y^3) 
Ω  = interpolate(ω, ProductSpace{GaussLobatto{M,18}, GaussLobatto{N,24}}) 
