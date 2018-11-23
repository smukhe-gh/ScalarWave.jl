#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
# Basis transformation Test
#--------------------------------------------------------------------

#---------------------------------------------
# test 1D basis transformation
#---------------------------------------------
ϕ = Field(GaussLobatto(U, 90), x->x^5 + exp(x) + sin(2x))
ψ = basistransform(ϕ)
λ = basistransform(ψ)
@test λ ≈ ϕ

#---------------------------------------------
# test 2D basis transformation
#---------------------------------------------

SUV = ProductSpace{GaussLobatto(U,2), GaussLobatto(V,4)}
𝕌 = Field(SUV, (U,V)->U)
𝕍 = Field(SUV, (U,V)->V)
𝕨 = exp(-((-5*𝕍^2 + 𝕌)^2)) 
𝕜 = exp(-((-5*𝕍^2 + 𝕌)^2)) 

# basis transformation using DFT 
𝕔_dft  = basistransform(𝕨, :dft) 
𝕨_dft  = basistransform(𝕔_dft, :dft)
@test 𝕨 ≈ 𝕨_dft

@assert 𝕨 == 𝕜
# basis transformation using MMT
𝕔_mmt = basistransform(𝕨)
𝕨_mmt = basistransform(𝕔_dft)
@test 𝕨 ≈ 𝕨_mmt
