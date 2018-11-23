#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 04-2018
# Test the derivative operators with ApproxFun
#--------------------------------------------------------------------

using ApproxFun, LinearAlgebra, HDF5

M  = 1.0
ω  = 1.0
l  = 0

PV, PU = 29, 29
Umax, Umin = -3M, -4M
Vmin, Vmax =  3M,  4M

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

# Test if the radial solution matches that of Mathematica
@test abs(ψ_re[50] - ψ_re_solved(expr[50])) < 1e-10

#--------------------------------------------------------------------
# Construct the 2D operator 
#--------------------------------------------------------------------
dU = -4M .. -3M; 
dV =  3M ..  4M;
d  = dU × dV
DU = ApproxFun.Derivative(d,[1,0]) 
DV = ApproxFun.Derivative(d,[0,1])

# Compute the r and t fields
r  = Fun((U,V) -> find_r_of_UV(U,V,1.0), d)
t  = Fun((U,V) -> find_t_of_UV(U,V,1.0), d)
invr = Fun((U,V) -> 1/find_r_of_UV(U,V,1.0), d)

# compute the operator
L  = DU*DV + ((DV*r)*invr)*DU + ((DU*r)*invr)*DV

#--------------------------------------------------------------------
# Set the boundary conditions 
#--------------------------------------------------------------------
U0 = Fun((_, V) -> ψ_re_solved(find_r_of_UV(-4M,  V, M)) * cos(ω * find_t_of_UV(-4M,  V, M)), ConstantSpace(ApproxFun.Point(-4M)) ⊗ ApproxFun.Chebyshev(3M .. 4M))
V0 = Fun((U, _) -> ψ_re_solved(find_r_of_UV(  U, 3M, M)) * cos(ω * find_t_of_UV(  U, 3M, M)), ApproxFun.Chebyshev(-4M .. -3M) ⊗ ConstantSpace(ApproxFun.Point(3M)))

B  = [I⊗ldirichlet(dV); ldirichlet(dU)⊗I]
u  = \([B; L], [V0; U0; 0]; tolerance=1E-12)

#--------------------------------------------------------------------
# Test with Mathematica 
#--------------------------------------------------------------------
ϕ_re = Fun((U,V)->ψ_re_solved(find_r_of_UV(U,V,M))*cos(ω * find_t_of_UV(U,V,M)), d)

@test_broken u(-3.3, 3.5) ≈ ϕ_re(-3.3, 3.5)
@test_broken u(-3.1, 3.9) ≈ ϕ_re(-3.1, 3.9)
@test_broken u(-3.4, 3.4) ≈ ϕ_re(-3.3, 3.4)

#--------------------------------------------------------------------
# Test with ScalarWave 
#--------------------------------------------------------------------

SUV = ScalarWave.ProductSpace{GaussLobatto(V,PV, Vmax, Vmin),
                   GaussLobatto(U,PU, Umax, Umin)}

𝕌 = Field(SUV, (U,V)->U)
𝕍 = Field(SUV, (U,V)->V)
t = Field(SUV, (U,V)->find_t_of_UV(U, V, M))
r = Field(SUV, (U,V)->find_r_of_UV(U, V, M))
ρ = 0 

𝔹 = boundary(Null, SUV)
𝔻𝕍, 𝔻𝕌 = derivative(SUV) 
𝔻r, 𝔻t = derivativetransform(SUV, t, r) 

ϕ_real = Field(SUV, (U,V) -> ψ_re_solved(find_r_of_UV(U,V,M)) * cos(ω * find_t_of_UV(U,V,M))) 
𝕓 = boundary(Null, SUV)*ϕ_real
𝕃 = 𝔻𝕌*𝔻𝕍 + ((𝔻𝕌*r)/r)*𝔻𝕍 +((𝔻𝕍*r)/r)*𝔻𝕌

# XXX: Scale the boundary operator
scaling = (1/(((PV+1)^2)*((PU+1)^2)))
invscaling = 1/scaling

@show cond(𝕃 + invscaling*𝔹)

# Compute the complex solution
𝕨 = solve(𝕃 + invscaling*𝔹, ρ + invscaling*𝕓) 

# Now compare the solutions; first the boundaries
u_collocation = Field(SUV, (U,V)->u(U,V))
testU =  𝕌.value[15, 24]
testV =  𝕍.value[15, 24]

@assert PV == PU
@test maximum(abs(ϕ_re(testU, testV) - ϕ_real.value[15, 24])) < 1e-10
for i in 1:PV
    @test maximum(abs(u(𝕌.value[i, 1], 𝕍.value[1, 1]) - ϕ_re(𝕌.value[i, 1], 𝕍.value[1, 1]))) < 1e-14 # boundary spanning U
    @test maximum(abs(u(𝕌.value[1, 1], 𝕍.value[1, i]) - ϕ_re(𝕌.value[1, 1], 𝕍.value[1, i]))) < 1e-14 # boundary spanning V
    @test maximum(abs(u(𝕌.value[i, 1], 𝕍.value[1, 1]) - ϕ_real.value[i, 1])) < 1e-14 # boundary spanning U
    @test maximum(abs(u(𝕌.value[1, 1], 𝕍.value[1, i]) - ϕ_real.value[1, i])) < 1e-14 # boundary spanning V
end

# compare the solutions
@show maximum(abs(u_collocation - 𝕨)) 
@show maximum(abs(ϕ_real - 𝕨)) 

# Compute the L2 norms
@show L2Error(𝕨, ϕ_real)
@show L2ErrorRelative(𝕨, ϕ_real)

#--------------------------------------------------------------------
# test if the operators satisfy the solution
#--------------------------------------------------------------------
@show maximum(abs(𝕃*ϕ_real))
