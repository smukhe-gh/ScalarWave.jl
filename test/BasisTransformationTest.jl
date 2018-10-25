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
𝕔_rfft = basistransform(𝕨)
𝕨_rfft = basistransform(𝕔_rfft)

# basis transformation using MMT
𝕨_mmt  = basistransform(𝕔_rfft, :MMT)
𝕔_mmt  = basistransform(𝕨_mmt,  :MMT) 

@test 𝕨 ≈ 𝕨_rfft
@test 𝕨 ≈ 𝕨_mmt
@test 𝕨_rfft ≈ 𝕨_mmt
@test 𝕔_rfft ≈ 𝕔_mmt

#---------------------------------------------
# Compare performance
#---------------------------------------------

using BenchmarkTools

for NN in 2:100
    SUV = ProductSpace{GaussLobatto{U,NN}, GaussLobatto{V,NN}}
    𝕌 = Field(SUV, (U,V)->U)
    𝕍 = Field(SUV, (U,V)->V)
    𝕨 = exp(-(𝕌^2 + 𝕍^2)) 

    if NN > 2
        println("DFT", @btime 𝕔_rfft = basistransform(𝕨))
        println("MMT", @btime 𝕔_mmt  = basistransform(𝕨_mmt,  :MMT))
    else
        𝕔_rfft = basistransform(𝕨)
        𝕔_mmt  = basistransform(𝕨_mmt,  :MMT) 
    end

end 
