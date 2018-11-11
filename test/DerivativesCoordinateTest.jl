#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 04-2018
# Test derivatives
#--------------------------------------------------------------------

using Einsum

#--------------------------------------------------------------------
# Define boundary and the product space
#--------------------------------------------------------------------
M = 1.0
ω = 1.0
PV, PU = 21, 21

Umin, Umax = -4M, -3M
Vmin, Vmax =  3M,  4M

SUV = ProductSpace{GaussLobatto(V,PV, Vmax, Vmin), GaussLobatto(U,PU, Umax, Umin)}

#--------------------------------------------------------------------
# Define coordinates and their associated derivatives
#--------------------------------------------------------------------
𝕌 = Field(SUV, (U,V)->U)
𝕍 = Field(SUV, (U,V)->V)

drawpatch(𝕌, "../output/scattering/coordinates/U")
drawpatch(𝕍, "../output/scattering/coordinates/V")

#--------------------------------------------------------------------
# Define derivative and boundary operators
#--------------------------------------------------------------------
𝔻𝕍, 𝔻𝕌 = derivative(SUV) 

#--------------------------------------------------------------------
# test the operator computation in isolation
# The issue seems to be with the coordinate transformation
#--------------------------------------------------------------------
r = sin(𝕌)*cos(𝕍)
𝕊 = (r^3)*(𝕌^4)*(𝕍^5)

#--------------------------------------------------------------------
# Now construct the operator according to 
# Carsten Gundlach and Jorge Pullin 1997 Class. Quantum Grav. 14 991
#--------------------------------------------------------------------
𝕃0  = 𝔻𝕌*𝔻𝕍 + ((𝔻𝕌*r)/r)*𝔻𝕍 +((𝔻𝕍*r)/r)*𝔻𝕌

𝕊0 = ( (𝕌^3)*(𝕍^4)*(cos(𝕍)^2)*(sin(𝕌)^2) * 
      (20*cos(𝕍)*(𝕌*cos(𝕌) + sin(𝕌)) - 𝕍*sin(𝕍)*(15*𝕌*cos(𝕌) + 16*sin(𝕌))) )

𝕃1 = 𝔻𝕌*𝔻𝕍
𝕃2 = 𝔻𝕌
𝕃3 = 𝔻𝕍

𝕊2 = (𝕌^3)*(𝕍^5)*(cos(𝕍)^3)*(sin(𝕌)^2)*(3*𝕌*cos(𝕌) + 4*sin(𝕌))
𝕊3 = (𝕌^4)*(𝕍^4)*(sin(𝕌)^3)*(cos(𝕍)^2)*(5*cos(𝕍) - 3*𝕍*sin(𝕍))
𝕊1 = (𝕌^3)*(𝕍^4)*(cos(𝕍)^2)*(sin(𝕌)^2)*(3*𝕌*cos(𝕌) + 4*sin(𝕌))*(5*cos(𝕍) - 3*𝕍*sin(𝕍))

@show maximum(abs(𝕃2*𝕊 - 𝕊2))
@show maximum(abs(𝕃3*𝕊 - 𝕊3))
@show maximum(abs(𝕃1*𝕊 - 𝕊1))
@show maximum(abs(𝕃0*𝕊 - 𝕊0))

# Hack to see how coefficents drop
ℂ = 1e-15 + basistransform(𝕊) 
@show minimum(abs(ℂ))
@show minimum(log10(abs(ℂ)))

using Plots
pyplot()
A = log10.(abs.(ℂ.value))
heatmap(A)
savefig("../output/2D-coefficents.pdf")
close()
