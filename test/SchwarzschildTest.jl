#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 04-2018
# Wave equation on Schwarzschild
#--------------------------------------------------------------------

using Einsum

struct U end
struct V end
struct UV end

tic()
#--------------------------------------------------------------------
# Define boundary and the product space
#--------------------------------------------------------------------
P1, P2 = 20, 20
M   = 1.0
Umin, Umax = -3M, -7M
Vmin, Vmax =  3M,  7M

#--------------------------------------------------------------------
# Define derivative and boundary operators
#--------------------------------------------------------------------
SUV = ProductSpace{GaussLobatto{U,P1}, GaussLobatto{V,P2}}
𝔹 = boundary(Null, SUV)

#--------------------------------------------------------------------
# Define coordinates and their associated derivatives
#--------------------------------------------------------------------
𝕌 = Field(SUV, (U,V)->U)
𝕍 = Field(SUV, (U,V)->V)
θ = Field(SUV, (𝑼,𝑽)->pi/2)
ϕ = Field(SUV, (𝑼,𝑽)->0)
ø = zero(SUV) 
Ø = zero(Spatial, SUV) 

𝑼 = (Umax + Umin)/2 + (Umax - Umin)/2*𝕌  
𝑽 = (Vmax + Vmin)/2 - (Vmax - Vmin)/2*𝕍  

t = Field(SUV, (𝑼,𝑽)->find_t_of_UV(𝑼, 𝑽, M), 𝑼, 𝑽)
r = Field(SUV, (𝑼,𝑽)->find_r_of_UV(𝑼, 𝑽, M), 𝑼, 𝑽)

𝔻𝑼, 𝔻𝑽 = derivativetransform(SUV, 𝑼, 𝑽) 
𝔻θ, 𝔻ϕ = Ø, Ø

#--------------------------------------------------------------------
# Set boundary conditions
#--------------------------------------------------------------------
ρ = 0 
𝕤 = exp(-((-5M + 𝑽)^2)) 
𝕓 = 𝔹*𝕤
drawpatch(𝕤, "boundary-field")

#--------------------------------------------------------------------
# Define metric functions 
#--------------------------------------------------------------------
𝒈𝑼𝑽 = -32*(M^3/r)*(exp(-r/2M))
𝒈θθ = r^2
𝒈ϕϕ = (r*sin(θ))^2

𝒈𝑼𝑼 = 𝒈𝑽𝑽 = ø
𝒈𝑼θ = 𝒈𝑼ϕ = ø
𝒈𝑽θ = 𝒈𝑽ϕ = ø
𝒈θϕ = ø

𝕘    = Metric{dd, 4}([𝒈𝑼𝑼, 𝒈𝑼𝑽, 𝒈𝑼θ, 𝒈𝑼ϕ, 
                           𝒈𝑽𝑽, 𝒈𝑽θ, 𝒈𝑽ϕ,
                                𝒈θθ, 𝒈θϕ,
                                     𝒈ϕϕ])

𝕘inv = metricinverse(𝕘) 
𝔻    = Derivative{u, 4}([𝔻𝑼, 𝔻𝑽, 𝔻θ, 𝔻ϕ])
Γ    = Christoffel(𝕘)
@einsum Γ[m, i, j] = (1/2)*𝕘inv[m,k]*(𝔻[j]*𝕘[k,i]+  𝔻[i]*𝕘[k,j] - 𝔻[k]*𝕘[i,j])
println("Finished computing auxilliary quantities. Computing the operator")
toc()
tic()
#--------------------------------------------------------------------
# Now construct the operator in 2 ways (just because you can)
#--------------------------------------------------------------------
𝕃1 = ( sum(𝕘inv[j,k]*𝔻[j]*𝔻[k] for j in 1:dim(𝕘), k in 1:dim(𝕘))  
     - sum(𝕘inv[j,k]*Γ[l,j,k]*𝔻[l] for j in 1:dim(𝕘), k in 1:dim(𝕘), l in 1:dim(𝕘)) ) 

#--------------------------------------------------------------------
# Solve the system [also check the condition number and eigen values]
#--------------------------------------------------------------------
𝕨 = solve(𝕃1 + 𝔹, ρ + 𝕓) 

println("Finished solve")
toc()

using PyPlot
plot(chebgrid(P2), 𝕨.value[:, end])
plot(chebgrid(P1), 𝕨.value[end, :])
plot(chebgrid(P1), 𝕨.value[1, :])
plot(chebgrid(P1), 𝕨.value[:, 1])
savefig("plots/boundary-R.pdf")

drawpatch(𝑼, "plots/transformedU")
drawpatch(𝑽, "plots/transformedV")
drawpatch(t, "plots/t_of_UV")
drawpatch(r, "plots/r_of_UV")
drawpatch(𝕨, "plots/schwarzschild")