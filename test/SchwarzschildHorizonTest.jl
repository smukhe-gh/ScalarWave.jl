#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 04-2018 Modified 01-2019
# Wave equation on Schwarzschild
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Define boundary and the product space
#--------------------------------------------------------------------
M = 1.0
ω = 1.0
PV, PU = 29, 29
Umax, Umin = (1/5)*M, -(4/5)*M
Vmin, Vmax =      3M,       4M
SUV = ProductSpace{GaussLobatto(V, PV, Vmax, Vmin),
                   GaussLobatto(U, PU, Umax, Umin)}

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
𝔻𝕍, 𝔻𝕌 = derivative(SUV) 

𝔻r, 𝔻t = derivativetransform(SUV, t, r) 
𝔻θ, 𝔻ϕ = Ø, Ø
ρ = 0 

#--------------------------------------------------------------------
# Set boundary conditions 
#--------------------------------------------------------------------
ϕ = exp((-(-(Vmax + Vmin)/2 + 𝕍)^2)/0.5) 
𝕓 = boundary(Null, SUV)*ϕ

#--------------------------------------------------------------------
# Now construct the operator according to 
# Carsten Gundlach and Jorge Pullin 1997 Class. Quantum Grav. 14 991
#--------------------------------------------------------------------
𝕃  = 𝔻𝕌*𝔻𝕍 + ((𝔻𝕌*r)/r)*𝔻𝕍 +((𝔻𝕍*r)/r)*𝔻𝕌

#--------------------------------------------------------------------
# Solve the system
#--------------------------------------------------------------------
𝕨 = solve(𝕃 + 𝔹, ρ + 𝕓) 
drawpatch(real(𝕨), "../output/real-psi")
