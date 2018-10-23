#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
# Basis transformation Test
#--------------------------------------------------------------------

ϕ = Field(GaussLobatto{U, 9}, x->x^5 + 2)
ψ = basistransform(ϕ, Chebyshev{U, 9})
λ = basistransform(ψ, GaussLobatto{U, 9})
