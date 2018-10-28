#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 04-2018
# Wave equation on Schwarzschild
#--------------------------------------------------------------------

using Einsum

#--------------------------------------------------------------------
# Define boundary and the product space
#--------------------------------------------------------------------
P1, P2 = 3, 3
M = 1.0
ω = 0.1 
Umin, Umax = -3M, -7M
Vmin, Vmax =  3M,  7M

#--------------------------------------------------------------------
# Define derivative and boundary operators
#--------------------------------------------------------------------
SUV = ProductSpace{GaussLobatto{U,P1}, GaussLobatto{V,P2}}
𝔹   = boundary(Null, SUV)

#--------------------------------------------------------------------
# Define coordinates and their associated derivatives
#--------------------------------------------------------------------
𝕌 = Field(SUV, (U,V)->U)
𝕍 = Field(SUV, (U,V)->V)
θ = Field(SUV, (𝑼,𝑽)->pi/2)
ϕ = Field(SUV, (𝑼,𝑽)->0)
ø = zero(SUV) 
Ø = zero(Null, SUV) 

𝑼 = (Umax + Umin)/2 + (Umax - Umin)/2*𝕌  
𝑽 = (Vmax + Vmin)/2 - (Vmax - Vmin)/2*𝕍  

t = Field(SUV, (𝑼,𝑽)->find_t_of_UV(𝑼, 𝑽, M), 𝑼, 𝑽)
r = Field(SUV, (𝑼,𝑽)->find_r_of_UV(𝑼, 𝑽, M), 𝑼, 𝑽)

𝔻𝑼, 𝔻𝑽 = derivativetransform(SUV, 𝑼, 𝑽) 
𝔻θ, 𝔻ϕ = Ø, Ø

#--------------------------------------------------------------------
# Set boundary conditions 
# [choose a solution and it's complex conjugate]
#--------------------------------------------------------------------
ρ = 0 

# Let the solution be of the form Exp(iω r)*Exp( -iω t)
# and it's conjugate (ω => -ω , t => -t)
# with unit amplitude, where '-iω' corresponds to
# the incoming wave

𝕤1re = cos(ω * r - ω * t)
s1im = sin(ω * r - ω * t)

𝕤2re = cos(-ω * r - ω * t)
s2im = sin(-ω * r - ω * t)

# Consider only an incoming wave
𝕓1re = boundary(Null, :R, SUV)*𝕤1re
𝕓1im = boundary(Null, :R, SUV)*s1im
𝕓2re = boundary(Null, :R, SUV)*𝕤2re
𝕓2im = boundary(Null, :R, SUV)*s2im

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

𝕘    = Metric{_dd, 4}([𝒈𝑼𝑼, 𝒈𝑼𝑽, 𝒈𝑼θ, 𝒈𝑼ϕ, 
                           𝒈𝑽𝑽, 𝒈𝑽θ, 𝒈𝑽ϕ,
                                𝒈θθ, 𝒈θϕ,
                                     𝒈ϕϕ])

𝕘inv = inv(𝕘) 
𝔻    = Derivative{_u, 4}([𝔻𝑼, 𝔻𝑽, 𝔻θ, 𝔻ϕ])
Γ    = Christoffel(𝕘)

@einsum Γ[m, i, j] = (1/2)*𝕘inv[m,k]*(𝔻[j]*𝕘[k,i]+  𝔻[i]*𝕘[k,j] - 𝔻[k]*𝕘[i,j])

#--------------------------------------------------------------------
# Now construct the operator in 2 ways (just because you can)
#--------------------------------------------------------------------
𝕃1 = ( sum(𝕘inv[j,k]*𝔻[j]*𝔻[k] for j in 1:dim(𝕘), k in 1:dim(𝕘))  
     - sum(𝕘inv[j,k]*Γ[l,j,k]*𝔻[l] for j in 1:dim(𝕘), k in 1:dim(𝕘), l in 1:dim(𝕘)) ) 

#--------------------------------------------------------------------
# Solve the system [also check the condition number and eigen values]
#--------------------------------------------------------------------
𝕨1re = solve(𝕃1 + 𝔹, ρ + 𝕓1re) 
𝕨1im = solve(𝕃1 + 𝔹, ρ + 𝕓1im) 

𝕨2re = solve(𝕃1 + 𝔹, ρ + 𝕓2re) 
𝕨2im = solve(𝕃1 + 𝔹, ρ + 𝕓2im) 

#--------------------------------------------------------------------
# Compute the Wronskian and check 
#--------------------------------------------------------------------

# Compute derivatives with respect to r
𝔻r, 𝔻t = derivativetransform(SUV, t, r)

dr_𝕨1re = 𝔻r*𝕨1re
dr_𝕨2re = 𝔻r*𝕨2re
dr_𝕨1im = 𝔻r*𝕨1im
dr_𝕨2im = 𝔻r*𝕨2im

# Wronskian for the real and the imaginary part; independently
Wre = 𝕨1re*dr_𝕨2re - 𝕨2re*dr_𝕨1re
Wim = 𝕨1im*dr_𝕨2im - 𝕨2im*dr_𝕨1im

# Now compute derivatives with respect to r to check if it's satisfied 
dr_Wre = 𝔻r*Wre
dr_Wim = 𝔻r*Wim

@show dr_Wre
@show dr_Wim
