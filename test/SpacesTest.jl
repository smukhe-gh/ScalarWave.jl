#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2018
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# 1D Spaces
#--------------------------------------------------------------------

struct M end

P = 5
S = Taylor{M, P}
ϕ = Field(S, x->x^3)  
γ = Field(S, x->x^5/20)  
b = Boundary(S, x->1//20, x->(-1//20))
D = derivative(S) 
B = boundary(S)
L = D*D
u = solve(L + B, ϕ + b)

@testset "1D spaces" begin
@test order(S) == P
@test dim(S)   == 1
@test range(S) == CartesianRange((P+1,))
@test size(S)  == (P+1,)
@test (ϕ + γ).value == Field(S, x->(x^5)/20 + x^3).value
@test (γ - ϕ).value == Field(S, x->(x^5)/20 - x^3).value
@test (γ * ϕ).value == Field(S, x->(x^5/20)*x^3).value
@test (D*γ).value   == Field(S, x->x^4/4).value
@test (ϕ*D*γ).value == Field(S, x->(x^3)*(x^4/4)).value
@test (L*γ).value   == Field(S, x->x^3).value
@test u.value       == γ.value 
end;

#--------------------------------------------------------------------
# 2D Spaces
#--------------------------------------------------------------------

struct U end
struct V end
struct UV end

SU  = GaussLobatto{U,3}
SV  = GaussLobatto{V,5}
SUV = ProductSpace{SU, SV}
ϕ   = Field(SUV, (x,y)->x+y)  
DU, DV = derivative(SUV)
B   = boundary(SUV)
L   = DV*DV + DU*DU


