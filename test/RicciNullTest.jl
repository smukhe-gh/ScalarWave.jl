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
M          = 1.0
P1, P2     = 20, 20
Umin, Umax = -3M, -7M
Vmin, Vmax =  3M,  7M

SUV = ProductSpace{GaussLobatto{U,P1}, GaussLobatto{V,P2}}

#--------------------------------------------------------------------
# Define coordinates and their associated derivatives
#--------------------------------------------------------------------
𝕌 = Field(SUV, (U,V)->U)
𝕍 = Field(SUV, (U,V)->V)
θ = Field(SUV, (𝑼,𝑽)->pi/2)
ϕ = Field(SUV, (𝑼,𝑽)->0)
ø = zero(SUV) 
Ø = zero(Spatial, SUV) 

𝑼 = (Umax + Umin)/2 + (Umax - Umin)/2*𝕌  
𝑽 = (Vmax + Vmin)/2 - (Vmax - Vmin)/2*𝕍  

t = Field(SUV, (𝑼,𝑽)->find_t_of_UV(𝑼, 𝑽, M), 𝑼, 𝑽)
r = Field(SUV, (𝑼,𝑽)->find_r_of_UV(𝑼, 𝑽, M), 𝑼, 𝑽)

𝔻𝑼, 𝔻𝑽 = derivativetransform(SUV, 𝑼, 𝑽) 
𝔻θ, 𝔻ϕ = Ø, Ø

#--------------------------------------------------------------------
# Define metric functions 
#--------------------------------------------------------------------
𝒈𝑼𝑽 = -32*(M^3/r)*(exp(-r/2M))
𝒈θθ = r^2
𝒈ϕϕ = (r*sin(θ))^2

𝒈𝑼𝑼 = 𝒈𝑽𝑽 = ø
𝒈𝑼θ = 𝒈𝑼ϕ = ø
𝒈𝑽θ = 𝒈𝑽ϕ = ø
𝒈θϕ = ø

𝕘   = Metric{dd, 4}([𝒈𝑼𝑼, 𝒈𝑼𝑽, 𝒈𝑼θ, 𝒈𝑼ϕ, 
                          𝒈𝑽𝑽, 𝒈𝑽θ, 𝒈𝑽ϕ,
                               𝒈θθ, 𝒈θϕ,
                                    𝒈ϕϕ])
𝕘inv = metricinverse(𝕘) 
𝔻    = [𝔻𝑼, 𝔻𝑽, 𝔻θ, 𝔻ϕ] 
Γ    = Christoffel(𝕘)
ℝ    = Ricci(𝕘)
@einsum Γ[m, i, j] = (1/2)*𝕘inv[m,k]*(𝔻[j]*𝕘[k,i]+  𝔻[i]*𝕘[k,j] - 𝔻[k]*𝕘[i,j])

#------------------------------------------------------
# Test Christoffels
#------------------------------------------------------

Γ111 = - (2M + r)*(𝔻𝑼*r)/(2M*r)
Γ222 = - (2M + r)*(𝔻𝑽*r)/(2M*r)
Γ133 = exp(r/2M)*(r^2)*(𝔻𝑽*r)/(32*M^3)
Γ441 = (𝔻𝑼*r)/r

@testset "Γ[a,b,c]" begin
    @test maximum(abs(Γ[1,1,1] - Γ111)) < 1e-10 
    @test maximum(abs(Γ[2,2,2] - Γ222)) < 1e-10 
    @test maximum(abs(Γ[1,3,3] - Γ133)) < 1e-10
    @test maximum(abs(Γ[4,4,1] - Γ441)) < 1e-10
end

quit()

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

