#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
# Basis transformation Test
#--------------------------------------------------------------------

#---------------------------------------------
# test 1D basis transformation
#---------------------------------------------
ϕ = Field(GaussLobatto(U, 9), x->x^5 + 2)
ψ = basistransform(ϕ)
λ = basistransform(ψ)

#---------------------------------------------
# test 2D basis transformation
#---------------------------------------------

Umin, Umax = -3, -7
Vmin, Vmax =  3,  7
SUV = ProductSpace{GaussLobatto(U,20), GaussLobatto(V,40)}

𝕌 = Field(SUV, (U,V)->U)
𝕍 = Field(SUV, (U,V)->V)
𝕨 = exp(-((-5*𝕍^2 + 𝕌)^2)) 

# basis transformation using MMT
𝕔_mmt = basistransform(𝕨)
𝕨_mmt = basistransform(𝕔_mmt)

# basis transformation using DFT 
𝕨_dft  = basistransform(𝕔_mmt, :dft)
𝕔_dft  = basistransform(𝕨_mmt,  :dft) 

drawpatch(𝕨, "w-field")
drawpatch(𝕨_mmt, "wmmt-field")
@test 𝕨_dft ≈ 𝕨_mmt
@test 𝕔_dft ≈ 𝕔_mmt
@test 𝕨 ≈ 𝕨_mmt
@test 𝕨 ≈ 𝕨_dft
exit()

#---------------------------------------------
# Test interpolation 
#---------------------------------------------

exit()
𝕏 = Field(ProductSpace{GaussLobatto(U,10), GaussLobatto(V,14)}, (U,V)->U + V) 
ℂ = basistransform(𝕏)
𝔻 = basistransform(ℂ)
@test 𝕏.value ≈ 𝔻.value

ℤ = interpolate(𝕏, ProductSpace{GaussLobatto(U,10), GaussLobatto(V,14)})
drawpatch(𝕏, "x-field")
drawpatch(ℤ, "z-field")

