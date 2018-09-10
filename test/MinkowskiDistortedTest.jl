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

#--------------------------------------------------------------------
# Define derivative and boundary operators
#--------------------------------------------------------------------
𝔻v, 𝔻u = derivative(SUV)
𝔹 = boundary(nullboundary, SUV)

testfield = Field(SUV,(u,v) -> u^2 + v^3)
@test (𝔻u*testfield).value ≈ Field(SUV,(u,v) -> 2*u).value
@test (𝔻v*testfield).value ≈ Field(SUV,(u,v) -> 3*v^2).value
testB = zeros(P2+1, P1+1)
testB[1, : ] = testB[:, 1] = 1.0
@test vec(𝔹) ≈ diagm(vec(testB))

#--------------------------------------------------------------------
# Set boundary conditions
#--------------------------------------------------------------------
𝕓 = Boundary(SUV, u->exp(-u^2/0.1), v->exp(-v^2/0.1), u->0, v->0)
ρ = Field(SUV, (u,v)->0) 

#--------------------------------------------------------------------
# Construct the wave operator in curved spacetime
#--------------------------------------------------------------------

# First we compute the auxiliary quantities
ω(u,v) = -(pi/8)*cospi(u/2)*cospi(v/2) 
Ω = Field(SUV, (u,v)->ω(u,v)) 
𝕤 = Field(SUV, (u,v)->exp(-(u*cos(ω(u,v)) + v*sin(ω(u,v)))^2/0.1) + exp(-(-u*sin(ω(u,v)) + v*cos(ω(u,v)))^2/0.1))

guu  = Field(SUV, (u,v)-> -4*cos(ω(u,v))*sin(ω(u,v)))
guv  = gvu =  Field(SUV, (u,v)-> -2*cos(2*ω(u,v)))
gvv  = Field(SUV, (u,v)-> 4*cos(ω(u,v))*sin(ω(u,v)))
detg = Field(SUV, (u,v)-> -1/4)

drawpatch(𝕤,    "plots/s-field")
#drawpatch(Ω,    "plots/omega-field")
#drawpatch(guu,  "plots/guu-field")
#drawpatch(guv,  "plots/guv-field")
#drawpatch(gvv,  "plots/gvv-field")

𝕃 = guu*𝔻u*𝔻u + gvv*𝔻v*𝔻v
  + (((2*detg)*(𝔻u*guv + 𝔻v*gvu) + (𝔻u*detg)*guu + (𝔻v*detg)*gvu)/(2*detg))*𝔻u
  + (((2*detg)*(𝔻u*guv + 𝔻v*gvv) + (𝔻u*detg)*guu + (𝔻v*detg)*gvv)/(2*detg))*𝔻v
  + (guu + gvu)*𝔻v*𝔻u

# A decent way (Erik's previous suggestion) to check the operator without a solve 
s = (𝕃 + 𝔹)*𝕤
@show maximum(s.value)
drawpatch(s, "plots/ss-field")

#--------------------------------------------------------------------
# Solve the system [also check the condition number and eigen values]
#--------------------------------------------------------------------
# TODO: Check various resolutions
#       Compare with analytic solution
#       Try with small twists first

# FIXME: We have a singular matrix with several zero eigen values. 

# XXX: Do we also need derivatives at the boundary?
𝔻𝔹 = (𝔻u + 𝔻v)*𝔹

# XXX: Setting Ω = 0, gives errors. This is most likely due to
# the wrong construction of the operator
@test_broken 𝕃 == -2*𝔻u*𝔻v;

# FIXME: We have zero eigen values and Inf condition number
@show sort(abs.(eigvals(𝕃 + 𝔹)))[1:10]
@show cond(𝕃 + 𝔹)

𝕨 = solve(𝕃 + 𝔹, ρ + 𝕓) 
drawpatch(𝕨, "plots/minkowski-distorted")
