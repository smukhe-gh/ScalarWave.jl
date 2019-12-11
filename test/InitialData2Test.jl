#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2019
# Simulate scalar field collpase on axis
#--------------------------------------------------------------------
# [5] Solve for an unknown spacetime. 

PS = ProductSpace(ChebyshevGL{U, 18, Float64}(-4, -2), 
                  ChebyshevGL{V, 18, Float64}( 4,  6))

#--------------------------------------------------------------------
# test with Schwarzschild spacetime first
#--------------------------------------------------------------------

M = 1
η = Field(PS, (u,v)->find_r_of_UV(u ,v, M))
f = (32*M^3/η)*exp(-η/2M)
a = sqrt(f)
ϕ = Field(PS, (u,v)->0)

(s, ψ, ϕ) = rescale(a, η, ϕ)
r = Field(PS, (u,v)->(v-u))

ronUbnd = extractUboundary(r, :incoming)
sonUbnd = extractUboundary(s, :incoming)
ψonUbnd = extractUboundary(ψ, :incoming)
ϕonUbnd = extractUboundary(ϕ, :incoming)
DV = derivative(sonUbnd.space)
I  = identity(sonUbnd.space)

@show L2(L_for_s(ψonUbnd, ϕonUbnd, ronUbnd, DV, I)*sonUbnd + rhs_for_s(ψonUbnd, ϕonUbnd, ronUbnd, DV))
sonUbndSolved = s_on_ubnd(s, ψ, ϕ)
@show L2(sonUbndSolved - sonUbnd)
@show L2(rescaledC2onbnd(sonUbndSolved, ψonUbnd, ϕonUbnd, ronUbnd, DV))

ψonUbndSolved = nonlinearsolver_for_psi(s, ψ, ϕ)
@show L2(ψonUbndSolved - ψonUbnd)

#--------------------------------------------------------------------
#Now try arbitrary spacetime
#--------------------------------------------------------------------

PS2 = ProductSpace(ChebyshevGL{U, 22, Float64}(0, 1), 
                   ChebyshevGL{V, 22, Float64}(0, 1))

r = Field(PS2, (u,v)->(v-u))
s = r
ψ = Field(PS2, (u,v)->1)
ϕ = Field(PS2, (u,v)->sin(u+v))

ronUbnd = extractUboundary(r, :incoming)
sonUbnd = extractUboundary(s, :incoming)
ψonUbnd = extractUboundary(ψ, :incoming)
ϕonUbnd = extractUboundary(ϕ, :incoming)
DV = derivative(sonUbnd.space)
I  = identity(sonUbnd.space)

sonUbndSolved = s_on_ubnd(s, ψ, ϕ)
plot(sonUbndSolved)
show()

# FIXME: Why doesn't this work on axis?
@show L2(rescaledC2onbnd(sonUbndSolved, ψonUbnd, ϕonUbnd, ronUbnd, DV))

# Let's try solving for ψ using the nonlinear solver
ψonUbndSolved = nonlinearsolver_for_psi(s, ψ, ϕ)
@show L2(rescaledC2onbnd(sonUbnd, ψonUbndSolved, ϕonUbnd, ronUbnd, DV))

