#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2018
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# 1D Spaces
#--------------------------------------------------------------------

struct M end
S = GaussLobatto{M, 9}
ϕ = Field(S, x->x+1)  
D = derivative(S) 

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
b   = Boundary(SUV, (x,y)->1)
DU, DV = derivative(SUV)
B   = boundary(SUV)
L   = DV*DV + DU*DU
