#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 06-2019
# Test datatypes and composability
#--------------------------------------------------------------------
struct M end
using LinearAlgebra

struct X end
struct Y end

SX  = ChebyshevGL{X, 3, BigFloat}(-1, 1)
SY  = ChebyshevGL{Y, 6, BigFloat}(-3, 5)
SXY = ProductSpace(SX, SY)

f = Field(SXY, (x,y)->x^2 +y^3)
g = Field(SXY, (x,y)->2x)
h = Field(SXY, (x,y)->3y^2)

DX = derivative(SX)
DY = derivative(SY)

D2X, D2Y = derivative(SXY)
@test (D2X*f).value ≈ g.value
@test (D2Y*f).value ≈ h.value


SW = ChebyshevGL{X, 10, BigFloat}(-1, pi)
DW = derivative(SW)
w  = Field(SW, x->x^2)
z  = Field(SW, x->2x)
@test DW*w ≈ z
