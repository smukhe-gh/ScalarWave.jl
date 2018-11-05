#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 04-2018
# Test covariant derivative 
#--------------------------------------------------------------------

using Einsum

#--------------------------------------------------------------------
# Define boundary and the product space
#--------------------------------------------------------------------
P1, P2 = 20, 20
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

𝒈tt = -(1 - 2M/𝒓)   
𝒈rr = 1/(1 - 2M/𝒓)  
𝒈θθ = 𝒓^2           
𝒈ϕϕ = (𝒓*sin(θ))^2  
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
R    = Ricci(𝕘)

@einsum Γ[m, i, j] = (1/2)*𝕘inv[m,k]*(𝔻[j]*𝕘[k,i]+  𝔻[i]*𝕘[k,j] - 𝔻[k]*𝕘[i,j])

#------------------------------------------------------
# Test Christoffels
#------------------------------------------------------

@testset "Γ[a,b,c]" begin
@test maximum(abs(Γ[1,1,2] - (M/𝒓^2)*((1 - 2(M/𝒓))^(-1)) )) < 1e-10  # 𝒕𝒕𝒓 
@test maximum(abs(Γ[2,1,1] - (M/𝒓^2)* (1 - 2(M/𝒓))       )) < 1e-10  # 𝒓𝒕𝒕
@test maximum(abs(Γ[2,2,2] + (M/𝒓^2)*((1 - 2(M/𝒓))^(-1)) )) < 1e-10  # 𝒓𝒓𝒓
@test maximum(abs(Γ[2,3,3] + (-2*M + 𝒓)                  )) < 1e-10  # 𝒓θθ
@test maximum(abs(Γ[2,4,4] + (-2*M + 𝒓)*(sin(θ)^2)       )) < 1e-10  # 𝒓ϕϕ
@test maximum(abs(Γ[3,2,3] - 1/𝒓                         )) < 1e-10  # θ𝒓θ
@test maximum(abs(Γ[4,2,4] - 1/𝒓                         )) < 1e-10  # ϕ𝒓ϕ
@test maximum(abs(Γ[3,4,4] + cos(θ)*sin(θ)               )) < 1e-10  # θϕϕ
@test maximum(abs(Γ[4,3,4] - (cos(θ)/sin(θ))             )) < 1e-10  # ϕθϕ

indices = ([1,1,2], [2,1,1], [2,2,2], [2,3,3],
           [2,4,4], [3,2,3], [4,2,4],
           [3,4,4], [4,3,4])

for a in 1:4, b in 1:4, c in 1:4
    @test Γ[a, b, c] == Γ[a, c, b]
    if !(([a,b,c] in indices) || ([a,c,b] in indices))
        @test maximum(abs(Γ[a,b,c])) < 1e-10
    end
end

end

#------------------------------------------------------
# Test Ricci 
#------------------------------------------------------

function computeRicci(𝔻, 𝕘, i, j)
    return (sum( 𝔻[l]*Γ[l,i,j] for l in 1:dim(𝕘) ) - 
            sum( 𝔻[j]*Γ[l,i,l] for l in 1:dim(𝕘) ) + 
            sum( Γ[m,i,j]*Γ[l,l,m] for m in 1:dim(𝕘), l in 1:dim(𝕘)) -  
            sum( Γ[m,i,l]*Γ[l,j,m] for m in 1:dim(𝕘), l in 1:dim(𝕘)) )
end

for i in 1:4, j in 1:4
    ℝ[i,j] = computeRicci(𝔻, 𝕘, i, j)
end

@testset "ℝ[a,b]" begin
    for i in 1:4, j in 1:4
        @test maximum(abs(ℝ[i,j])) < 1e-8
    end
end

# NOTE: We should have two broken tests. 
#       For R[3,3] and R[4,4] since we do not take the derivatives correctly

#------------------------------------------------------
# Construct covariant derivatives and compatibility 
#------------------------------------------------------



















