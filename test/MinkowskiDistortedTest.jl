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
#--------------------------------------------------------------------
nullboundary = Null
P1, P2 = 4, 4
SUV = ProductSpace{GaussLobatto{U,P1}, GaussLobatto{V,P2}}

#@show boundary(y, SUV).value

#--------------------------------------------------------------------
# Define derivative and boundary operators
#--------------------------------------------------------------------
𝔻v, 𝔻u = derivative(SUV)
𝔹 = boundary(nullboundary, SUV)

#--------------------------------------------------------------------
# Set boundary conditions
#--------------------------------------------------------------------
𝕓 = Boundary(SUV, u->exp(-u^2/0.1), v->0*exp(-v^2/0.1), u->0, v->0)
ρ = Field(SUV, (u,v)->0) 

#--------------------------------------------------------------------
# Construct the wave operator in curved spacetime
#--------------------------------------------------------------------

# First we compute the auxiliary quantities
ω(u,v) = (pi/8)*cospi(u/2)*cospi(v/2) 
Ω = Field(SUV, (u,v)->ω(u,v)) 
𝕤 = Field(SUV, (u,v)->exp(-(u*cos(ω(u,v)) + v*sin(ω(u,v)))^2/0.1) + exp(-(-u*sin(ω(u,v)) + v*cos(ω(u,v)))^2/0.1))

guu  = Field(SUV, (u,v)-> -4*cos(ω(u,v))*sin(ω(u,v)))
guv  = gvu =  Field(SUV, (u,v)-> -2*cos(2*ω(u,v)))
gvv  = Field(SUV, (u,v)-> 4*cos(ω(u,v))*sin(ω(u,v)))
detg = Field(SUV, (u,v)-> -1/4)

𝕘 = [guu guv; gvu gvv]
𝔻 = [𝔻u, 𝔻v]
det𝕘 = detg

Ł = sum((𝕘[a,b]*𝔻[a])*𝔻[b] + ((sqrt(1/abs(det𝕘))*𝔻[a]*(𝕘[a,b]*sqrt(abs(det𝕘)))))*𝔻[b] for a in 1:2, b in 1:2)
𝕃 = (guu*𝔻u*𝔻u + gvv*𝔻v*𝔻v + guv*𝔻u*𝔻v + gvu*𝔻v*𝔻u
     + sqrt(1/abs(detg))*(𝔻u*(guu*sqrt(abs(detg))) + 𝔻v*(gvu*sqrt(abs(detg))))*𝔻u
     + sqrt(1/abs(detg))*(𝔻u*(guv*sqrt(abs(detg))) + 𝔻v*(gvv*sqrt(abs(detg))))*𝔻v)

@test maximum(𝕃.value - Ł.value) < 1e-13

#drawpatch(𝕤,    "plots/s-field")
#drawpatch(Ω,    "plots/omega-field")
#drawpatch(guu,  "plots/guu-field")
#drawpatch(guv,  "plots/guv-field")
#drawpatch(gvv,  "plots/gvv-field")

#--------------------------------------------------------------------
# Solve the system [also check the condition number and eigen values]
#--------------------------------------------------------------------

𝕨 = solve(𝕃 + 𝔹, ρ + 𝕓) 
drawpatch(𝕨, "plots/minkowski-distorted")
