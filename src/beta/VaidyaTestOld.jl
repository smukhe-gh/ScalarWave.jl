#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 12-2018
# Represent the analytic Vaidya metric and do the tests 
#   -- Compute the expansion at each point
#   -- Compute the Weyl scalar ψ_2 
#--------------------------------------------------------------------

using Roots, Einsum

# define the constants for an exponential mass function
α = 1.0
β = 1.0
c = 1.0

# construct the functions r[u,v] and f[u,v]
function m(v)
    return (1/β)*(α*exp((β/2)*c*v) + 1)
end

function P(u)
    return -c*u
end

function x(r, v)
    return sqrt(r^2 - (4r/β) + (4*m(v)/β))
end

function r(u,v)
    f(r) = β*x(r, v) + 2*log(r - (2/β) + x(r,v)) - (β*c*v)/2 - P(u)
    r    = find_zero(f, 2)
    return r
end

function f(u,v)
    return (-c^2*x(u,v))/(β*r(u,v))
end

# now construct the space and the operator
PV, PU = 29, 29
Vmax, Vmin = 10, 2
Umax, Umin = -30, -60
SUV = ScalarWave.ProductSpace{GaussLobatto(V,PV, Vmax, Vmin),
                              GaussLobatto(U,PU, Umax, Umin)}
𝔻𝕍, 𝔻𝕌 = derivative(SUV) 
𝑟 = Field(SUV, (u,v)->r(u,v))
𝑓 = Field(SUV, (u,v)->f(u,v))
θ = Field(SUV, (u,v)->π/2)
ø = zero(SUV) 

drawpatch(𝑟, "../output/r-coordinate-vaidya")
drawpatch(𝑓, "../output/f-coordinate-vaidya")

# Now compute the expansion at each point
# See Poisson for the details of the derivation
Θ = 2*(𝔻𝕍*𝑟)/(𝑟*𝑓)

drawpatch(Θ, "../output/expansion-vaidya")

# test the metric quantities
g = Metric{_dd, 4}([ø, 𝑓, ø, ø, 
                       ø, ø, ø,
                          𝑟, ø,
                            (𝑟^2)*sin(θ)^2])

invg = 
𝔻    = Derivative{_u, 4}([𝔻𝕌, 𝔻𝕍, 𝔻θ, 𝔻ϕ])
Γ    = Christoffel(g)

# compute Chirstoffels
@einsum Γ[m, i, j] = (1/2)*invg[m,k]*(𝔻[j]*g[k,i]+  𝔻[i]*g[k,j] - 𝔻[k]*g[i,j])

# compute Ricci
function computeRicci(𝔻, 𝕘, i, j)
    return (sum( 𝔻[l]*Γ[l,i,j] for l in 1:dim(𝕘) ) - 
            sum( 𝔻[j]*Γ[l,i,l] for l in 1:dim(𝕘) ) + 
            sum( Γ[m,i,j]*Γ[l,l,m] for m in 1:dim(𝕘), l in 1:dim(𝕘)) -  
            sum( Γ[m,i,l]*Γ[l,j,m] for m in 1:dim(𝕘), l in 1:dim(𝕘)) )
end

# test components of the Ricci tensor with Mathematica/Waugh and Lake
R11 = computeRicci(𝔻, 𝗀, 1, 1)
R22 = computeRicci(𝔻, 𝗀, 2, 2)
R33 = computeRicci(𝔻, 𝗀, 3, 3)
R44 = computeRicci(𝔻, 𝗀, 4, 4)
R14 = computeRicci(𝔻, 𝗀, 1, 4)


