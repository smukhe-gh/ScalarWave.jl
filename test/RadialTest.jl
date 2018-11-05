#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 04-2018
# Test the radial wave operator with Mathematica
# --------------------------------------------------------------------

using HDF5

#--------------------------------------------------------------------
# Construct the radial operator and test with Mathematica 
# interpolating function
#--------------------------------------------------------------------
M  = 1.0
ω  = 1.0
l  = 0
#--------------------------------------------------------------------
# Choose and intermediate domain close to the black hole 
#--------------------------------------------------------------------
rmax, rmin = 200M, 3M 
SR = GaussLobatto(U, 1000, rmax, rmin) 
𝔻r = derivative(SR) 
I  = eye(SR)
r  = Field(SR, r->r)
f  = 1 - (2M/r)

@time 𝕃  = (f^2)*𝔻r*𝔻r + (2M/r^2)*f*𝔻r + (ω^2 - f*( (2M/r^3) + (l*(l+1)/(r^2)) ))*I


# Export data to Mathematica to compute the interpolation functioin at the 
if isfile("../output/collocation-points-for-mathematica.h5")
    println("File already exists.")
else
    h5write("../output/collocation-points-for-mathematica.h5", "collocation-points", r.value)
end

# Import the data from Mathematica and load it into an array
if isfile("../output/values-for-julia.h5")
    ψ_re = Field(SR, h5read("../output/values-for-julia.h5", "psi-real"))
    ψ_im = Field(SR, h5read("../output/values-for-julia.h5", "psi-imaginary"))
else
    println("Waiting for Mathematica to generate files")
    exit()
end

residual_re = 𝕃*ψ_re
residual_im = 𝕃*ψ_im

using Plots
pyplot()
plot(r.value[10:990],  residual_re.value[10:990], leg=false)
plot!(r.value[10:990], residual_im.value[10:990], leg=false)
savefig("../output/residual-plot-$rmin-$rmax.pdf")

plot(r.value[10:90],  residual_re.value[10:90], leg=false)
plot!(r.value[10:90], residual_im.value[10:90], leg=false)
savefig("../output/residual-plot-zoomed-$rmin-$rmax.pdf")

plot(r.value[10:990],  ψ_re.value[10:990], leg=false)
plot!(r.value[10:990], ψ_im.value[10:990], leg=false)
savefig("../output/psi-plot-$rmin-$rmax.pdf")
