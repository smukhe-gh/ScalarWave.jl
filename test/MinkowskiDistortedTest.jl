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
P1, P2 = 40, 40
SUV = ProductSpace{GaussLobatto{U,P1}, GaussLobatto{V,P2}}

#--------------------------------------------------------------------
# Define derivative and boundary operators
#--------------------------------------------------------------------
𝔻v, 𝔻u = derivative(SUV)
𝔹 = boundary(nullboundary, SUV)

#--------------------------------------------------------------------
# Define coordinates and their associated derivatives
#--------------------------------------------------------------------
u = Field(SUV, (u,v)->u)
v = Field(SUV, (u,v)->v)
Ω = Field(SUV, (u,v)->(pi/8)*cospi(u/2)*cospi(v/2))

𝒖 =  u*cos(Ω) + v*sin(Ω)
𝒗 = -u*sin(Ω) + v*cos(Ω)

𝔻𝒖ofu = ( cos(-Ω) + (-sin(-Ω)*u + cos(-Ω)*v)*(𝔻u*(-Ω))) 
𝔻𝒖ofv = ( sin(-Ω) + (-sin(-Ω)*u + cos(-Ω)*v)*(𝔻v*(-Ω))) 
𝔻𝒗ofu = (-sin(-Ω) + (-cos(-Ω)*u - sin(-Ω)*v)*(𝔻u*(-Ω))) 
𝔻𝒗ofv = ( cos(-Ω) + (-cos(-Ω)*u - sin(-Ω)*v)*(𝔻v*(-Ω))) 

𝔻𝒖    = 𝔻𝒖ofu * 𝔻u + 𝔻𝒖ofv * 𝔻v  
𝔻𝒗    = 𝔻𝒗ofu * 𝔻u + 𝔻𝒗ofv * 𝔻v

#========================================
 Test the chain rule for derivatives
========================================#

G = sin(𝒖)
H = 𝒖^2

@test (𝔻u*𝒖).value ≈ ( cos(Ω) + (-sin(Ω)*u + cos(Ω)*v)*(𝔻u*Ω)).value 
@test (𝔻v*𝒖).value ≈ ( sin(Ω) + (-sin(Ω)*u + cos(Ω)*v)*(𝔻v*Ω)).value 
@test (𝔻u*𝒗).value ≈ (-sin(Ω) + (-cos(Ω)*u - sin(Ω)*v)*(𝔻u*Ω)).value 
@test (𝔻v*𝒗).value ≈ ( cos(Ω) + (-cos(Ω)*u - sin(Ω)*v)*(𝔻v*Ω)).value 

@test (𝔻u*G).value ≈ (cos(𝒖)*(cos(Ω) + (-sin(Ω)*u + cos(Ω)*v)*(𝔻u*Ω))).value 
@test (𝔻u*H).value ≈ ((2*𝒖)*(cos(Ω) + (-sin(Ω)*u + cos(Ω)*v)*(𝔻u*Ω))).value 

# SubTest: Test if the Jacobian computation is accurate
𝔻uof𝒖 = 𝔻u*𝒖
𝔻vof𝒖 = 𝔻v*𝒖
𝔻uof𝒗 = 𝔻u*𝒗
𝔻vof𝒗 = 𝔻v*𝒗

𝔻𝒖ofu = Field(SUV, similar(𝔻uof𝒖.value)) 
𝔻𝒖ofv = Field(SUV, similar(𝔻vof𝒖.value)) 

𝔻𝒗ofu = Field(SUV, similar(𝔻uof𝒗.value))
𝔻𝒗ofv = Field(SUV, similar(𝔻vof𝒗.value))

# Check if reading in and reading out the Jacobian is correct.
for index in CartesianRange(size(𝔻uof𝒖.value)) 
    Jacobian = [𝔻uof𝒖.value[index] 𝔻uof𝒗.value[index]; 
                𝔻vof𝒖.value[index] 𝔻vof𝒗.value[index]]
    InverseJacobian    = inv(Jacobian)
    𝔻𝒖ofu.value[index] = InverseJacobian[1,1] 
    𝔻𝒖ofv.value[index] = InverseJacobian[1,2] 
    𝔻𝒗ofu.value[index] = InverseJacobian[2,1] 
    𝔻𝒗ofv.value[index] = InverseJacobian[2,2] 
end

𝔻𝒖    = 𝔻𝒖ofu * 𝔻u + 𝔻𝒖ofv * 𝔻v  
𝔻𝒗    = 𝔻𝒗ofu * 𝔻u + 𝔻𝒗ofv * 𝔻v

@test (𝔻𝒖*G).value ≈ cos(𝒖).value 

quit()

#--------------------------------------------------------------------
# Set boundary conditions
#--------------------------------------------------------------------

ρ = Field(SUV, (u,v)->0)
𝕤 = exp(-(𝒖^2)/0.1) 
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

