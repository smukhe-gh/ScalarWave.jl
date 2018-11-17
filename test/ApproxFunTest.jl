#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 04-2018
# Test the derivative operators with ApproxFun
#--------------------------------------------------------------------

using ApproxFun, LinearAlgebra, HDF5

M  = 1.0
ω  = 1.0
l  = 0

#=
#--------------------------------------------------------------------
# Solve the radial equation 
#--------------------------------------------------------------------
dR = 3M .. 20M
r  = Fun(r->r, dR)
invr = Fun(r->1/r, dR)
f  = 1 - (2M*invr)
Dr = ApproxFun.Derivative(dR)

# compute operator and boundary conditioins
L  = f^2*Dr*Dr + (2M*(invr^2)*f)*Dr + (ω^2 - f*( (2M*(invr^3)) + (l*(l+1)*(invr^2)) ))*ApproxFun.I
B  = Dirichlet(dR)

# Import the data from Mathematica and load it into an array
if isfile("../output/hdf5/values-for-julia.h5")
    ψ_re = h5read("../output/hdf5/values-for-julia.h5", "psi-real")
    ψ_im = h5read("../output/hdf5/values-for-julia.h5", "psi-imag")
    expr = h5read("../output/hdf5/collocation-points-for-mathematica.h5", "collocation-points")
else
    println("Waiting for Mathematica to generate files")
    exit()
end

ψ_re_solved = [B;L] \ [[ψ_re[1], ψ_re[end]], 0]
@test abs(ψ_re[50] - ψ_re_solved(expr[50])) < 1e-10
=#

#--------------------------------------------------------------------
# Compute the 2D operator 
#--------------------------------------------------------------------

dU = -4M .. -3M; 
dV =  3M ..  4M;
d  = dU × dV
DU = ApproxFun.Derivative(d,[1,0]) 
DV = ApproxFun.Derivative(d,[0,1])

#=
# Compute the r and t fields
r  = Fun((U,V) -> find_r_of_UV(U,V,1.0), d)
t  = Fun((U,V) -> find_t_of_UV(U,V,1.0), d)
invr = Fun((U,V) -> 1/find_r_of_UV(U,V,1.0), d)

L  = DU*DV + ((DV*r)*invr)*DU + ((DU*r)*invr)*DV

#--------------------------------------------------------------------
# Set the boundary conditions 
#--------------------------------------------------------------------

# Compute the solution
ϕ_re = Fun((U,V)->ψ_re_solved(find_r_of_UV(U,V,M))*cos(ω * find_t_of_UV(U,V,M)), d)

U_bnd = Fun(V->ψ_re_solved(find_r_of_UV(-4M, 1.0*V, M)) * cos(ω * find_t_of_UV(-4M, 1.0*V, M)), dV)
V_bnd = Fun(U->ψ_re_solved(find_r_of_UV(1.0*U, 3M, M))  * cos(ω * find_t_of_UV(1.0*U, 3M, M)),  dU)

SUV = ScalarWave.ProductSpace{GaussLobatto(V, 150, 4M, 3M), GaussLobatto(U, 150, -3M, -4M)}
𝕌   = Field(SUV, (U,V)->U)
𝕍   = Field(SUV, (U,V)->V)

# for testing if the boundary conditions are applied correctly
ϕ_re_array = zeros(151, 151)
for _u in 1:151, _v in 1:151
    ϕ_re_array[_u, _v] = ϕ_re(𝕌.value[_u, _v], 𝕍.value[_u, _v]) 
end

# Now check if the boundaries are being computed properly
"""
using Plots
pyplot()
plot(U_bnd)
plot!(V_bnd)
savefig("../output/boundary-UV.pdf")
close()
"""
=#

# Test
V0 = Fun(x->exp(-20x^2),dU)
L  = DU*DV

# TODO: scale the boundary operator appropriately
B  = [I⊗ldirichlet(dV); ldirichlet(dU)⊗I]
u  = \([B; L], [V0; 0; 0]; tolerance=1E-12)

#=
@test u(-3.3, 3.5) ≈ ϕ_re(-3.3, 3.5)
@test u(-3.1, 3.9) ≈ ϕ_re(-3.1, 3.9)
@test u(-3.4, 3.4) ≈ ϕ_re(-3.3, 3.4)
println("Finished Testing solution")
=#
