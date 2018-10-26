#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
# Basis transformation Test
#--------------------------------------------------------------------

#---------------------------------------------
# test 1D basis transformation
#---------------------------------------------
ϕ = Field(GaussLobatto{U, 9}, x->x^5 + 2)
ψ = basistransform(ϕ, Chebyshev{U, 9})
λ = basistransform(ψ, GaussLobatto{U, 9})

#---------------------------------------------
# test 2D basis transformation
#---------------------------------------------

Umin, Umax = -3, -7
Vmin, Vmax =  3,  7
SUV = ProductSpace{GaussLobatto{U,2}, GaussLobatto{V,4}}

𝕌 = Field(SUV, (U,V)->U)
𝕍 = Field(SUV, (U,V)->V)
𝑼 = (Umax + Umin)/2 + (Umax - Umin)/2*𝕌  
𝑽 = (Vmax + Vmin)/2 - (Vmax - Vmin)/2*𝕍  

𝕨 = exp(-((-5 + 𝑽)^2)) 

# basis transformation using DFT
𝕔_mmt = basistransform(𝕨)
𝕨_mmt = basistransform(𝕔_mmt)

# basis transformation using MMT
𝕨_dft  = basistransform(𝕔_mmt, :dft)
𝕔_dft  = basistransform(𝕨_mmt,  :dft) 

@test 𝕨 ≈ 𝕨_mmt
@test 𝕨 ≈ 𝕨_dft
@test 𝕨_dft ≈ 𝕨_mmt
@test 𝕔_dft ≈ 𝕔_mmt
