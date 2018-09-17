#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2018
# Distorted Minkowski
#--------------------------------------------------------------------

struct U end
struct V end
struct UV end
struct 𝑼𝑽 end 

#--------------------------------------------------------------------
# Define boundary and the product space
# Derivative tests fails for P <= 20
#--------------------------------------------------------------------
nullboundary = Null
P1, P2 = 40, 40
SUV = ProductSpace{GaussLobatto{U,P1}, GaussLobatto{V,P2}}
S𝑼𝑽 = ProductSpace{GaussLobatto{U,P1}, GaussLobatto{V,P2}}

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

𝔻𝒗, 𝔻𝒖 = derivative(SUV, 𝒖, 𝒗)

G = sin(𝒖)
H = 𝒖^2
@test (𝔻𝒖*G).value ≈ cos(𝒖).value 
@test (𝔻𝒖*𝒖).value ≈ Field(SUV, (u,v)->1).value 
@test maximum((𝔻𝒖*𝒗).value) < 1e-12
@test_broken (𝔻𝒖*𝒗).value ≈ Field(SUV, (u,v)->0).value 

#--------------------------------------------------------------------
# Set boundary conditions
#--------------------------------------------------------------------

ρ = 0*sin(𝒖 + 𝒗)*exp(-𝒖^2/0.1)*exp(-𝒗^2/0.1)
𝕤 = exp(-(𝒗^2)/0.1) 
𝕓 = 𝔹*𝕤

#--------------------------------------------------------------------
# Construct the wave operator in curved spacetime
#--------------------------------------------------------------------

𝕘uu  = -4*cos(Ω)*sin(Ω)
𝕘uv  = 𝕘vu = -2*cos(2*Ω)
𝕘vv  = 4*cos(Ω)*sin(Ω)
det𝕘 = Field(SUV, (u,v)-> -1/4)

𝕘    = [𝕘uu 𝕘uv; 𝕘vu 𝕘vv]
𝔻    = [𝔻𝒖, 𝔻𝒗]
𝕃    = sum(( 𝕘[a,b] * 𝔻[a]) * 𝔻[b] 
           + ((sqrt(1/abs(det𝕘)) * 𝔻[a] * (𝕘[a,b]*sqrt(abs(det𝕘))))) * 𝔻[b] for a in 1:2, b in 1:2)

#--------------------------------------------------------------------
# Solve the system [also check the condition number and eigen values]
#--------------------------------------------------------------------
𝕨 = solve(𝕃 + 𝔹, ρ + 𝕓) 
drawpatch(𝕨, "plots/minkowski-distorted")

