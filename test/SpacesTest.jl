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
@test len(S)   == P+1
@test (ϕ + γ)  == Field(S, x->(x^5)/20 + x^3)
@test (γ - ϕ)  == Field(S, x->(x^5)/20 - x^3)
@test (γ * ϕ)  == Field(S, x->(x^5/20)*x^3)
@test (D*γ)    == Field(S, x->x^4/4)
@test (ϕ*D*γ)  == Field(S, x->(x^3)*(x^4/4))
@test (L*γ)    == Field(S, x->x^3)
@test u        == γ
end

P = 20
S = GaussLobatto{M, P}
ϕ = Field(S, x-> exp(4x))  
w = Field(S, x->(exp(4x) - x*sinh(4.0) - cosh(4.0))/16)  
b = Boundary(S, x->0, x->0)
D = derivative(S) 
B = boundary(S)
L = D*D
v = solve(L + B, ϕ + b)
@test v ≈ w 

#--------------------------------------------------------------------
# 2D Spaces
#--------------------------------------------------------------------

struct U end
struct V end
struct UV end

SS = Spatial
NN = Null

P1, P2 = 3, 5
SU  = Taylor{U,P1}
SV  = Taylor{V,P2}
SUV = ProductSpace{SU, SV}

γ   = Field(SUV, (x,y)->0)  
ϕ   = Field(SUV, (x,y)->x+y)
ψ   = Field(SUV, (x,y)->x^2+y^3)  
dψ  = Field(SUV, (x,y)->x^2+y^3+5)  
DU, DV = derivative(SUV)
B   = boundary(Spatial, SUV)
I   = identity(SU) ⦼ identity(SV)
b   = Boundary(SUV, x->x^2 + 1, y->1+y^3, x->x^2 - 1, y->1+y^3)

@testset "2D spaces" begin
@test order(SUV) == (P2, P1)
@test dim(SUV)   == 2
@test range(SUV) == CartesianRange((P2+1,P1+1))
@test size(SUV)  == (P2+1,P1+1)
@test (ϕ + ψ)  == Field(SUV, (x,y)->x^2 + y^3 + x + y) 
@test (ϕ - ψ)  == Field(SUV, (x,y)->x + y - x^2 - y^3) 
@test (ϕ * ψ)  == Field(SUV, (x,y)->(x^2 + y^3)*(x + y)) 
@test (ψ / dψ) == Field(SUV, (x,y)->(x^2 + y^3)/(x^2 + y^3 + 5))
@test DU.value == reshape(kron(derivative(SU).value, identity(SV).value), (P2+1,P1+1,P2+1,P1+1))
@test DV.value == reshape(kron(identity(SU).value, derivative(SV).value), (P2+1,P1+1,P2+1,P1+1))
@test (DU*DV).value == reshape(kron(derivative(SU).value, identity(SV).value)*
                               kron(identity(SU).value, derivative(SV).value), (P2+1,P1+1,P2+1,P1+1))
@test (ϕ*I).value == reshape(diagm(vec(ϕ.value)), (P2+1, P1+1, P2+1, P1+1))
@test (I*ψ) == ψ
@test (DU*ψ) == Field(SUV, (x,y)->3y^2)
@test (DU*DU*ψ) == Field(SUV, (x,y)->6y)
@test (DV*ψ) == Field(SUV, (x,y)->2x)
@test (DU*DV + DV*DU).value == (DV*DU + DU*DV).value
end;

P1, P2 = 5, 7
SU  = Taylor{U,P1}
SV  = Taylor{V,P2}
SUV = ProductSpace{SU, SV}
ψ   = Field(SUV, (x,y)->x^2+y^3 + x^3*y^2)  
dxψ = Field(SUV, (x,y)->2x + 3*x^2*y^2)  
dyψ = Field(SUV, (x,y)->3y^2 + 2*x^3*y)  
ddxddyψ = Field(SUV, (x,y)->2 + 2x^3 + 6y + 6x*y^2)

𝔹 = boundary(Spatial, SUV)
B = zeros(Rational{BigInt}, size(SUV))
B[1, :] = B[:, 1] = B[:, end] = B[end, :] = 1//1
b = Boundary(SUV, x->x^2 + x^3 + 1, y->y^3 + y^2 + 1, x->x^2 - 1 + x^3, y->1 + y^3 - y^2)

Dy, Dx = derivative(SUV)
Ł = Dx*Dx + Dy*Dy
𝕓 = 𝔹*ψ
𝕦 = reshape((Ł + 𝔹).value, (prod(size(SUV)), prod(size(SUV)))) \ vec((ddxddyψ + 𝕓).value)
𝕨 = solve(Ł + 𝔹, ddxddyψ + 𝕓) 

@testset "2D Laplace Solve" begin
@test (Dx*ψ) == dxψ 
@test (Dy*ψ) == dyψ 
@test reshape(dxψ.value, prod(size(SUV))) == vec(dxψ.value)
@test reshape(Dx.value, prod(size(SUV)), prod(size(SUV)))*vec(ψ.value) == vec(dxψ.value)
@test reshape(𝔹.value, prod(size(SUV)), prod(size(SUV))) == diagm(vec(B))
@test b.value == (𝔹*ψ).value

@test (Ł*ψ).value == ddxddyψ.value
@test 𝕦 == vec(ψ.value)
@test 𝕨 == ψ
end
