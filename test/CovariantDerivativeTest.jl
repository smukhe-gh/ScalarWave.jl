#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 04-2018
# Test Ricci and Weyl tensor computations on Schwarzschild
#--------------------------------------------------------------------

using Einsum

#--------------------------------------------------------------------
# Define boundary and the product space
#--------------------------------------------------------------------
P1, P2 = 30, 30
M = 1.0
SUV = ProductSpace{GaussLobatto(U,P1, 5M, 3M), 
                   GaussLobatto(V,P2, 5M, 3M)}

#--------------------------------------------------------------------
# Define coordinates and their associated derivatives
#--------------------------------------------------------------------

𝒕 = Field(SUV, (u,v)->u)
𝒓 = Field(SUV, (u,v)->v)
θ = Field(SUV, (u,v)->pi/2)
ϕ = Field(SUV, (u,v)->0)

ø = zero(SUV) 
Ø = zero(Null, SUV) 

𝔻𝒓, 𝔻𝒕 = derivative(SUV) 
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

𝕘    = Metric{_dd, 4}([𝒈tt, 𝒈tr, 𝒈tθ, 𝒈tϕ, 
                           𝒈rr, 𝒈rθ, 𝒈tϕ,
                                𝒈θθ, 𝒈θϕ,
                                     𝒈ϕϕ])

𝕘inv = inv(𝕘) 
𝔻    = Derivative{_u, 4}([𝔻𝒕, 𝔻𝒓, 𝔻θ, 𝔻ϕ])

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
# Compute covariant derivatives and check metric
# compatibility
# FIXME: Why are these tests failing? 
#------------------------------------------------------

CD = Array{Union{Nothing, Field}}(nothing, 4,4,4) 

for a in 1:dim(𝕘), b in 1:dim(𝕘), c in 1:dim(𝕘)
    CD[a,b,c] = ( 𝔻[c]*𝕘[a,b]
                 + sum( Γ[m,c,a]*𝕘[m,b] for m in 1:dim(𝕘))
                 + sum( Γ[n,c,b]*𝕘[a,n] for n in 1:dim(𝕘)) )
end

@testset "CD[c]*g[a,b]" begin
for a in 1:dim(𝕘), b in 1:dim(𝕘), c in 1:dim(𝕘)
    @test maximum(abs(CD[a,b,c])) < 1e-8
end
end

