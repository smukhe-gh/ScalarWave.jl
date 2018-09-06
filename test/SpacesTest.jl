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

@testset "1D space" begin
@test order(S) == P
@test dim(S)   == 1
@test range(S) == 1:P+1
@test len(S)  == P+1
@test (ϕ + γ).value == Field(S, x->(x^5)/20 + x^3).value
@test (γ - ϕ).value == Field(S, x->(x^5)/20 - x^3).value
@test (γ * ϕ).value == Field(S, x->(x^5/20)*x^3).value
@test (D*γ).value   == Field(S, x->x^4/4).value
@test (ϕ*D*γ).value == Field(S, x->(x^3)*(x^4/4)).value
@test (L*γ).value   == Field(S, x->x^3).value
@test u.value       == γ.value 
end;

P = 20
S = GaussLobatto{M, P}
ϕ = Field(S, x->exp(4x))  
w = Field(S, x->(exp(4x) - x*sinh(4.0) - cosh(4.0))/16)  
b = Boundary(S, x->0, x->0)
D = derivative(S) 
B = boundary(S)
L = D*D
v = solve(L + B, ϕ + b)
@test v.value ≈ w.value 

#--------------------------------------------------------------------
# 2D Spaces
#--------------------------------------------------------------------


struct U end
struct V end
struct UV end

P1, P2 = 3, 5
SU  = Taylor{U,P1}
SV  = Taylor{V,P2}
SUV = ProductSpace{SU, SV}

ϕ   = Field(SUV, (x,y)->x+y)  
γ   = Field(SUV, (x,y)->0)  
ψ   = Field(SUV, (x,y)->x^2+y^3)  
DU, DV = derivative(SUV)
B   = boundary(SUV)
L   = DU*DU + DV*DV
I   = identity(SU) ⦼ identity(SV)
b   = Boundary(SUV, x->0//1, x->0//1, x->0//1, x->x^2) 

@testset "2D spaces" begin
@test order(SUV) == (P2, P1)
@test dim(SUV)   == 2
@test range(SUV) == CartesianRange((P2+1,P1+1))
@test size(SUV)  == (P2+1,P1+1)
@test (ϕ + ψ).value == Field(SUV, (x,y)->x^2 + y^3 + x + y).value 
@test (ϕ - ψ).value == Field(SUV, (x,y)->x + y - x^2 - y^3).value 
@test (ϕ * ψ).value == Field(SUV, (x,y)->(x^2 + y^3)*(x + y)).value 
@test DU.value == reshape(kron(derivative(SU).value, identity(SV).value), (P2+1,P1+1,P2+1,P1+1))
@test DV.value == reshape(kron(identity(SU).value, derivative(SV).value), (P2+1,P1+1,P2+1,P1+1))
@test (DU*DV).value == reshape(kron(derivative(SU).value, identity(SV).value)*
                               kron(identity(SU).value, derivative(SV).value), (P2+1,P1+1,P2+1,P1+1))
@test (ϕ*I).value == reshape(diagm(vec(ϕ.value)), (P2+1, P1+1, P2+1, P1+1))
@test (I*ψ).value == ψ.value 
@test (DU*ψ).value == Field(SUV, (x,y)->3y^2).value
@test (DU*DU*ψ).value == Field(SUV, (x,y)->6y).value
@test (DV*ψ).value == Field(SUV, (x,y)->2x).value
end;
