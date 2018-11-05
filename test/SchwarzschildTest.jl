#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 04-2018
# Wave equation on Schwarzschild
# Tests to be done
#   -- Test if the metric is compatible.
#   -- Compute the Weyl tensor to check for the fall-off conditions?
#   -- Use a different computation of the operator and check if
#      the operator constructions agree. 
#   -- Use the solution from Mathematica and check if our operator
#      satisfies the solution [Done]
#   -- Divide out the time-dependence (i.e. convert into a 
#      stationary solution) and check if the solution is independent 
#      of time. This is not a-priori obvious. 
#   -- Is is possible to plug this solution into the differential 
#      operator and check? We'd need knowledge of l and m.   
#--------------------------------------------------------------------

using Einsum

#--------------------------------------------------------------------
# Define boundary and the product space
#--------------------------------------------------------------------
M = 1.0
ω = 1.0
PV, PU = 40, 40
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
𝒓 = r + 2M*log(-1 + r/2M)

#--------------------------------------------------------------------
# Define derivative and boundary operators
#--------------------------------------------------------------------
𝔹 = boundary(Null, SUV)
𝔻𝕌, 𝔻𝕍 = derivative(SUV) 
𝔻r, 𝔻t = derivativetransform(SUV, t, r) 
𝔻θ, 𝔻ϕ = Ø, Ø

#--------------------------------------------------------------------
# Set boundary conditions 
# Note that you'd need to start with a set of boundary conditions
# that satisfy the operator.
#--------------------------------------------------------------------

ρ = 0 
𝕤 = exp(im*𝒓)*exp(-im*ω*t) 
𝕓 = boundary(SUV, Null, :R)*𝕤

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

#--------------------------------------------------------------------
# Now construct the operator 
#--------------------------------------------------------------------

𝕃1 = ( sum(𝕘inv[j,k]*𝔻[j]*𝔻[k] for j in 1:dim(𝕘), k in 1:dim(𝕘))  
     - sum(𝕘inv[j,k]*Γ[l,j,k]*𝔻[l] for j in 1:dim(𝕘), k in 1:dim(𝕘), l in 1:dim(𝕘)) ) 

#--------------------------------------------------------------------
# Solve the system [also check the condition number and eigen values]
#--------------------------------------------------------------------
𝕨 = solve(𝕃1 + 𝔹, ρ + 𝕓) 

#--------------------------------------------------------------------
# [T1] Check for time-stationarity 
#--------------------------------------------------------------------


