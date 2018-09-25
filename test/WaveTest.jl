#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2018
# Distorted Minkowski
#--------------------------------------------------------------------

struct U end
struct V end
struct UV end

#--------------------------------------------------------------------
# Define boundary and the product space
# Derivative tests fails for P <= 20
#--------------------------------------------------------------------
nullboundary = Null
P1, P2 = 40, 40
SUV = ProductSpace{GaussLobatto{U,P1}, GaussLobatto{V,P2}}

#--------------------------------------------------------------------
# Define derivative and boundary operators
#--------------------------------------------------------------------
𝔹 = boundary(nullboundary, SUV)

#--------------------------------------------------------------------
# Define coordinates and their associated derivatives
#--------------------------------------------------------------------
u = Field(SUV, (u,v)->u)
v = Field(SUV, (u,v)->v)
Ω = Field(SUV, (u,v)->(pi/8)*cospi(u/2)*cospi(v/2))

𝒖 =  u*cos(Ω) + v*sin(Ω)
𝒗 = -u*sin(Ω) + v*cos(Ω)
𝔻𝒗, 𝔻𝒖 = derivativetransform(SUV, 𝒖, 𝒗)

#--------------------------------------------------------------------
# Set boundary conditions
#--------------------------------------------------------------------
ρ = 0 
𝕤 = exp(-((𝒖^2)/0.1)) 
𝕓 = 𝔹*𝕤

#--------------------------------------------------------------------
# Construct the wave operator in curved spacetime
#--------------------------------------------------------------------
guu = Field(SUV, (u,v)-> 0)
guv = Field(SUV, (u,v)->-2)
gvv = Field(SUV, (u,v)-> 0)

(𝕘𝒖𝒖, 𝕘𝒖𝒗, 𝕘𝒗𝒗) = inversemetrictransform(guu, guv, gvv, 𝒖, 𝒗) 
invsqrtdet𝕘     = 1/sqrt(abs(inversemetricdet(𝕘𝒖𝒖, 𝕘𝒖𝒗, 𝕘𝒗𝒗))) 

𝕘   = [𝕘𝒖𝒖 𝕘𝒖𝒗; 𝕘𝒖𝒗 𝕘𝒗𝒗]
𝔻   = [𝔻𝒖, 𝔻𝒗]
𝕃   = 𝕘𝒖𝒗*𝔻𝒖*𝔻𝒗 + 𝕘𝒖𝒗*𝔻𝒗*𝔻𝒖

#--------------------------------------------------------------------
# Solve the system [also check the condition number and eigen values]
#--------------------------------------------------------------------
𝕨 = solve(𝕃 + 𝔹, ρ + 𝕓) 
drawpatch(𝕨, "plots/minkowski-distorted")
@show maximum(abs(𝕃*𝕤))
