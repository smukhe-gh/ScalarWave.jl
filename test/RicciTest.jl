#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 04-2018
# Test Ricci and Weyl tensor computations on Schwarzschild
#--------------------------------------------------------------------

struct U end
struct V end
struct UV end

#--------------------------------------------------------------------
# Define boundary and the product space
#--------------------------------------------------------------------
P1, P2 = 40, 40
M = 1.0
SUV = ProductSpace{GaussLobatto{U,P1}, GaussLobatto{V,P2}}

#--------------------------------------------------------------------
# Define coordinates and their associated derivatives
#--------------------------------------------------------------------

t = Field(SUV, (u,v)->u)
r = Field(SUV, (u,v)->v)
θ = Field(SUV, (u,v)->pi/2)
ϕ = Field(SUV, (u,v)->0)

ø = Field(SUV, (u,v)->0)
Ø = zero(SUV) 

𝒕 = ((4M - 3M)/2)*t  + (4M + 3M)/2
𝒓 = ((4M - 3M)/2)*r  + (4M + 3M)/2

𝔻𝒓, 𝔻𝒕 = derivativetransform(SUV, 𝒕, 𝒓) 
𝔻θ, 𝔻ϕ = Ø, Ø

#--------------------------------------------------------------------
# Define metric functions 
#--------------------------------------------------------------------

𝒈tt = (1 - 2M/𝒓)
𝒈rr = 1/(1 - 2M/𝒓)
𝒈θθ = 𝒓^2 
𝒈ϕϕ = (𝒓*sin(θ))^2 
𝒈rθ = 𝒈rϕ = 𝒈tr = ø 
𝒈tθ = 𝒈tϕ = 𝒈θϕ = ø

𝕘   = Metric{dd, 2}([𝒈tt; 
                  𝒈tr 𝒈rr; 
                  𝒈tθ 𝒈rθ 𝒈θθ;
                  𝒈tϕ 𝒈rϕ 𝒈θϕ 𝒈ϕϕ])

𝕘inv = metricinverse(𝕘) 

𝔻    = [𝔻𝒕 𝔻𝒓 𝔻θ 𝔻ϕ]

@einsum Γddd[k,i,j] := (1/2)*(𝔻[j]*𝕘[k,i] + 𝔻[i]*𝕘[k,j] - 𝔻[k]*𝕘[i,j])
@einsum Γ[m,i,j]    := 𝕘inv[m,k]*Γddd[k,i,j]



