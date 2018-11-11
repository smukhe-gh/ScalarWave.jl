#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 04-2018
# Wave equation on Schwarzschild
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Define boundary and the product space
#--------------------------------------------------------------------
M = 1.0
ω = 2.0
PV, PU = 29, 29
Umax, Umin = -3M, -3M
Vmin, Vmax =  3M,  3M
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
𝔻𝕍, 𝔻𝕌 = derivative(SUV) 

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
𝕃  = 𝔻𝕌*𝔻𝕍 + ((𝔻𝕌*r)/r)*𝔻𝕍 +((𝔻𝕍*r)/r)*𝔻𝕌

#--------------------------------------------------------------------
# Solve the system
#--------------------------------------------------------------------
𝕨 = solve(𝕃 + 𝔹, ρ + 𝕓) 
𝕔 = basistransform(real(𝕨)) + im*basistransform(imag(𝕨))

#--------------------------------------------------------------------
# Visualize solutions 
#--------------------------------------------------------------------

drawpatch(𝕌, "../output/scattering/coordinates/U")
drawpatch(𝕍, "../output/scattering/coordinates/V")
drawpatch(t, "../output/scattering/coordinates/t")
drawpatch(r, "../output/scattering/coordinates/r")
drawpatch(real(ϕ), "../output/scattering/waves/phi-real")
drawpatch(imag(ϕ), "../output/scattering/waves/phi-imag")
drawpatch(real(𝕨), "../output/scattering/waves/wave-real")
drawpatch(imag(𝕨), "../output/scattering/waves/wave-imag")

"""
using Plots
pyplot()
A = log10(abs(real(𝕔))).value
B = log10(abs(imag(𝕔))).value
heatmap(A)
savefig("../output/scattering/coeffs/coeffs_real.pdf")
heatmap(B)
savefig("../output/scattering/coeffs/coeffs_imag.pdf")
close()
"""

#--------------------------------------------------------------------
# Compare solutions 
#--------------------------------------------------------------------

@show maximum(abs(real(ϕ) - real(𝕨)))
@show maximum(abs(imag(ϕ) - imag(𝕨)))

drawpatch(abs(real(ϕ) - real(𝕨)), "../output/scattering/error/error-wave_real")
drawpatch(abs(imag(ϕ) - imag(𝕨)), "../output/scattering/error/error-wave_imag")

#--------------------------------------------------------------------
# Check for time-stationarity 
#--------------------------------------------------------------------
𝕧 = 𝕨 * exp(im * ω  * t) 
ψ = ϕ * exp(im * ω  * t) 
 
Dt_𝕧 = 𝔻t * 𝕧
Dt_ψ = 𝔻t * ψ

# Take time derivatives and check
@show maximum(abs(real(𝔻t * 𝕧)))
@show maximum(abs(imag(𝔻t * 𝕧)))

drawpatch(real(𝔻t * 𝕧), "../output/scattering/error/error-time-derivative-wave-real")
drawpatch(imag(𝔻t * 𝕧), "../output/scattering/error/error-time-derivative-wave-imag")
drawpatch(real(𝔻t * ψ), "../output/scattering/error/error-time-derivative-phi-real")
drawpatch(imag(𝔻t * ψ), "../output/scattering/error/error-time-derivative-phi-imag")

#--------------------------------------------------------------------
# Check what's happening at the boundaries 
#--------------------------------------------------------------------

u_bndOL = 𝕌.value[:, end]
v_bndOR = 𝕍.value[end, :]

ϕ_bndOL = ϕ.value[:, end]
𝕨_bndOL = 𝕨.value[:, end]
ϕ_bndOR = ϕ.value[end, :]
𝕨_bndOR = 𝕨.value[end, :]


using Plots
pyplot()
plot( u_bndOL, real(ϕ_bndOL), lab="phi-outgoing-left")
plot!(u_bndOL, real(𝕨_bndOL), lab="sol-outgoing-left",  line=:dot)
savefig("../output/scattering/boundaries/boundaries-u-real-outgoing-left.pdf")
close()

plot( v_bndOR, real(ϕ_bndOR), lab="phi-outgoing-right")
plot!(v_bndOR, real(𝕨_bndOR), lab="sol-outgoing-right", line=:dot)
savefig("../output/scattering/boundaries/boundaries-v-real-outgoing-right.pdf")
close()
