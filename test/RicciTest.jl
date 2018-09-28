#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 04-2018
# Test Ricci and Weyl tensor computations on Schwarzschild
#--------------------------------------------------------------------

using Einsum

struct U end
struct V end
struct UV end

#--------------------------------------------------------------------
# Define boundary and the product space
#--------------------------------------------------------------------
P1, P2 = 15, 15
M = 1.0
SUV = ProductSpace{GaussLobatto{U,P1}, GaussLobatto{V,P2}}

#--------------------------------------------------------------------
# Define coordinates and their associated derivatives
#--------------------------------------------------------------------

t = Field(SUV, (u,v)->u)
r = Field(SUV, (u,v)->v)
θ = Field(SUV, (u,v)->pi/2)
ϕ = Field(SUV, (u,v)->0)

ø = zero(SUV) 
Ø = zero(Spatial, SUV) 

𝒕 = (5M + 3M)/2 + ((5M - 3M)/2)*t  
𝒓 = (5M + 3M)/2 + ((5M - 3M)/2)*r  

𝔻𝒓, 𝔻𝒕 = derivativetransform(SUV, 𝒕, 𝒓) 
𝔻θ, 𝔻ϕ = Ø, Ø

#--------------------------------------------------------------------
# Define metric functions 
#--------------------------------------------------------------------

𝒈tt = -(1 - 2M/𝒓)    # Field(SUV, (u,v)->1) 
𝒈rr = 1/(1 - 2M/𝒓)  # Field(SUV, (u,v)->2) 
𝒈θθ = 𝒓^2           # Field(SUV, (u,v)->3) 
𝒈ϕϕ = (𝒓*sin(θ))^2  # Field(SUV, (u,v)->4) 
𝒈rθ = 𝒈rϕ = 𝒈tr = ø 
𝒈tθ = 𝒈tϕ = 𝒈θϕ = ø

𝕘    = Metric{dd, 4}([𝒈tt, 𝒈tr, 𝒈tθ, 𝒈tϕ, 
                           𝒈rr, 𝒈rθ, 𝒈tϕ,
                                𝒈θθ, 𝒈θϕ,
                                     𝒈ϕϕ])


𝕘inv = metricinverse(𝕘) 
𝔻    = Derivative{u, 4}([𝔻𝒕, 𝔻𝒓, 𝔻θ, 𝔻ϕ])

Γ    = Christoffel(𝕘)
ℝ    = Ricci(𝕘)
δ    = eye(4)

@einsum Γ[m, i, j] = (1/2)*𝕘inv[m,k]*(𝔻[j]*𝕘[k,i]+  𝔻[i]*𝕘[k,j] - 𝔻[k]*𝕘[i,j])
@einsum ℝ[i,j] = 𝔻[l]*Γ[l,i,j] - 𝔻[j]*Γ[l,i,l] + Γ[m,i,j]*Γ[l,l,m] - Γ[m,i,l]*Γ[l,j,m]
@einsum ℝ[i,j] = 𝔻[j]*Γ[l,i,l]


@testi ℝ[1,1] ==  (𝔻[1]*Γ[1,1,1] +
                   𝔻[1]*Γ[2,1,2] + 
                   𝔻[1]*Γ[3,1,3] +
                   𝔻[1]*Γ[4,1,4])


@test maximum(abs(Γ[2,1,1] - (M/𝒓^3)*(-2M + 𝒓))) < 1e-13
@test Γ[3,2,3] ≈ 1/𝒓
@test Γ[4,2,4] ≈ 1/𝒓

for i in 1:4, j in 1:4
    @show i, j, maximum(abs(ℝ[i,j]))
end


