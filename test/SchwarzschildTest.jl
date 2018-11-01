#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 04-2018
# Wave equation on Schwarzschild
#--------------------------------------------------------------------

using Einsum

#--------------------------------------------------------------------
# Define boundary and the product space
#--------------------------------------------------------------------
M = 1.0
PV, PU = 30, 30
Umax, Umin = -3M, -7M
Vmin, Vmax =  3M,  7M
SUV = ProductSpace{GaussLobatto(V,PV, Vmax, Vmin), 
                   GaussLobatto(U,PU, Umax, Umin)}

#--------------------------------------------------------------------
# Define coordinates and their associated derivatives
#--------------------------------------------------------------------
𝕌 = Field(SUV, (U,V)->U)
𝕍 = Field(SUV, (U,V)->V)
θ = Field(SUV, (U,V)->π/2)
ϕ = Field(SUV, (U,V)->0)

ø = zero(SUV) 
Ø = zero(Null, SUV) 

t = Field(SUV, (U,V)->find_t_of_UV(U, V, M), 𝕌, 𝕍)
r = Field(SUV, (U,V)->find_r_of_UV(U, V, M), 𝕌, 𝕍)

#--------------------------------------------------------------------
# Define derivative and boundary operators
#--------------------------------------------------------------------
𝔹 = boundary(Null, SUV)
𝔻𝕌, 𝔻𝕍 = derivative(SUV) 
𝔻θ, 𝔻ϕ = Ø, Ø

#--------------------------------------------------------------------
# Set boundary conditions
#--------------------------------------------------------------------
ρ = 0 
𝕤 = exp(-((-5M + 𝕍)^2)) 
𝕓 = 𝔹*𝕤

#--------------------------------------------------------------------
# Define metric functions 
#--------------------------------------------------------------------
𝒈𝕌𝕍 = -32*(M^3/r)*(exp(-r/2M))
𝒈θθ = r^2
𝒈ϕϕ = (r*sin(θ))^2

𝒈𝕌𝕌 = 𝒈𝕍𝕍 = ø
𝒈𝕌θ = 𝒈𝕌ϕ = ø
𝒈𝕍θ = 𝒈𝕍ϕ = ø
𝒈θϕ = ø

𝕘   = Metric{_dd, 4}([𝒈𝕌𝕌, 𝒈𝕌𝕍, 𝒈𝕌θ, 𝒈𝕌ϕ, 
                           𝒈𝕍𝕍, 𝒈𝕍θ, 𝒈𝕍ϕ,
                                𝒈θθ, 𝒈θϕ,
                                     𝒈ϕϕ])

𝕘inv = inv(𝕘) 
𝔻    = Derivative{_u, 4}([𝔻𝕌, 𝔻𝕍, 𝔻θ, 𝔻ϕ])
Γ    = Christoffel(𝕘)
@einsum Γ[m, i, j] = (1/2)*𝕘inv[m,k]*(𝔻[j]*𝕘[k,i]+  𝔻[i]*𝕘[k,j] - 𝔻[k]*𝕘[i,j])

# Comptute the Weyl tensor
# Compute the solution that's stationary in time
# and test the operator
#--------------------------------------------------------------------
# Now construct the operator 
# L1 is in terms of partial derivatives. L2 and L3 are both in terms
# of covariant derivatives, and at the very least, both should be
# equivalent since the first covariant derivative is a scalar.
#--------------------------------------------------------------------

𝕃1 = ( sum(𝕘inv[j,k]*𝔻[j]*𝔻[k] for j in 1:dim(𝕘), k in 1:dim(𝕘))  
     - sum(𝕘inv[j,k]*Γ[l,j,k]*𝔻[l] for j in 1:dim(𝕘), k in 1:dim(𝕘), l in 1:dim(𝕘)) ) 

#𝕃2 = sum( 𝕘[a,b]*(𝔻[b]*𝔻[a] - sum( Γ[k,a,b]*𝔻[k] for k in 1:dim(𝕘) )) for a in 1:dim(𝕘), b in 1:dim(𝕘) )
#𝕃3 = sum( 𝔻[b]*𝔻[b] + sum( Γ[b,b,k]*𝔻[k] for k in 1:dim(𝕘) ) for  b in 1:dim(𝕘) )

#@show maximum(abs.((𝕃1 - 𝕃2).value)) 
#@show maximum(abs(𝕃1 - 𝕃3)) 
#@show maximum(abs(𝕃2 - 𝕃3)) 

#--------------------------------------------------------------------
# Solve the system [also check the condition number and eigen values]
#--------------------------------------------------------------------
𝕨 = solve(𝕃1 + 𝔹, ρ + 𝕓) 
𝕔 = basistransform(𝕨)

@show maximum(abs(𝕃1*𝕨))
