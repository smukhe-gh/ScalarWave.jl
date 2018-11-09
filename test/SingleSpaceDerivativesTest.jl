#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 04-2018
# Test derivatives
#--------------------------------------------------------------------

PU   =  80
Umin = -4
Umax =  4
SU   = GaussLobatto(U,PU, Umax, Umin)

#--------------------------------------------------------------------
# Define coordinates and their associated derivatives
#--------------------------------------------------------------------
𝕌  = Field(SU, U->U)
𝔻𝕌 = derivative(SU)  

#--------------------------------------------------------------------
# test the operator computation in isolation
# The issue seems to be with the coordinate transformation
#--------------------------------------------------------------------
𝕊 = (sin(𝕌)^3)*(𝕌^4)

𝕃2 = 𝔻𝕌
𝕊2 = (𝕌^3)*(sin(𝕌)^2)*(3*𝕌*cos(𝕌) + 4*sin(𝕌))

coeffs = basistransform(𝕊)
@show maximum(abs(𝕃2*𝕊 - 𝕊2))


using Plots
pyplot()
plot(log10.(abs.(coeffs.value)), line=:dot)
savefig("coeffs.pdf")
close()

"""
using Plots
pyplot()
plot(𝕌.value, 𝕊.value)
plot!(𝕌.value, 𝕊2.value)
plot!(𝕌.value, (𝕃2*𝕊).value, line=:dot)
savefig("plotUS.pdf")
"""
