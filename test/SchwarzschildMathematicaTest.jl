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
ω = 1.0
PV, PU = 10, 10
Umax, Umin = -4M, -8M
Vmin, Vmax =  4M,  8M
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
ρ = 0 

#--------------------------------------------------------------------
# Set boundary conditions 
#--------------------------------------------------------------------
using HDF5

# Choose to export the schwarzschild r
if isfile("../output/hdf5/collocation-points-r.h5")
    println("File already exits. Skipping")
else
    println("Creating dataset.")
    h5open("../output/hdf5/collocation-points-r.h5", "w") do file
        write(file, "collocation-points-grid",  r.value)
    end
end

# Read ϕ(r) from Mathematica and multiply by exp(-iω t)
if isfile("../output/hdf5/values-for-julia-grid.h5")
    ϕr_real = Field(SUV, h5read("../output/hdf5/values-for-julia-grid.h5", "psi-real-grid"))
    ϕr_imag = Field(SUV, h5read("../output/hdf5/values-for-julia-grid.h5", "psi-imag-grid"))
else
    println("Could not find file. Create them using Mathematica")
    exit()
end

ϕ = (ϕr_real + im*ϕr_imag)*exp(-im * ω * t)
𝕓 = boundary(Null, SUV)*ϕ

#--------------------------------------------------------------------
# Now construct the operator according to 
# Carsten Gundlach and Jorge Pullin 1997 Class. Quantum Grav. 14 991
#--------------------------------------------------------------------
𝕃 = 𝔻𝕌*𝔻𝕍 + ((𝔻𝕌*r)/r)*𝔻𝕍 + ((𝔻𝕍*r)/r)*𝔻𝕌

#--------------------------------------------------------------------
# Solve the system [also check the condition number and eigen values]
#--------------------------------------------------------------------
𝕨 = solve(𝕃 + 𝔹, ρ + 𝕓) 

# FIXME: compute the coefficents: work with complex fields
# 𝕔 = basistransform(𝕨_real)

#--------------------------------------------------------------------
# Visualize solutions 
#--------------------------------------------------------------------
drawpatch(𝕌, "../output/scattering/U")
drawpatch(𝕍, "../output/scattering/V")
drawpatch(t, "../output/scattering/t")
drawpatch(r, "../output/scattering/r")

drawpatch(real(ϕ), "../output/scattering/phi-r-real")
drawpatch(imag(ϕ), "../output/scattering/phi-r-imag")
drawpatch(real(𝕨), "../output/scattering/wave_real")
drawpatch(imag(𝕨), "../output/scattering/wave_imag")

exit()

"""
using Plots
pyplot()
A = log(abs(𝕔_real)).value
heatmap(A)
savefig("../output/scattering/coefficents_real.pdf")
close()

B = log(abs(𝕔_imag)).value
heatmap(B)
savefig("../output/scattering/coefficents_imaginary.pdf")
close()
"""

#--------------------------------------------------------------------
# Compare solutions 
#--------------------------------------------------------------------

@show maximum(abs(ϕ_real - 𝕨_real))
@show maximum(abs(ϕ_imag - 𝕨_imag))

drawpatch(ϕ_real - 𝕨_real, "../output/scattering/error-wave_real")
drawpatch(ϕ_imag - 𝕨_imag, "../output/scattering/error-wave_imag")

#--------------------------------------------------------------------
# Check for time-stationarity 
#--------------------------------------------------------------------
𝕨_real_real = 𝕨_real*cos(ω * t)
𝕨_real_imag = 𝕨_real*sin(ω * t)

𝕨_imag_real = 𝕨_imag*cos(ω * t)
𝕨_imag_imag = 𝕨_imag*sin(ω * t)

𝕨_stationary_real = 𝕨_real_real - 𝕨_imag_imag
𝕨_stationary_imag = 𝕨_real_imag + 𝕨_imag_real 

# Take time derivatives and check
@show maximum(abs(𝔻t*𝕨_stationary_real))
@show maximum(abs(𝔻t*𝕨_stationary_imag))

