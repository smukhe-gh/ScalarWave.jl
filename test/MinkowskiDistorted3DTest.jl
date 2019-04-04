#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2018
# Distorted Minkowski
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Define boundary and the product space
# Derivative tests fails for P <= 20
#--------------------------------------------------------------------
nullboundary = Null
P1, P2, P3 = 2, 2, 2
SUVW = ProductSpace{GaussLobatto(U,P1), 
                    GaussLobatto(V,P2),
                    GaussLobatto(W,P3)}

#--------------------------------------------------------------------
# Define derivative and boundary operators
#--------------------------------------------------------------------
B = boundary(nullboundary, SUVW)

#--------------------------------------------------------------------
# Define coordinates and their associated derivatives
#--------------------------------------------------------------------
u = Field(SUVW, (u,v,w)->u)
v = Field(SUVW, (u,v,w)->v)
w = Field(SUVW, (u,v,w)->w)

Ω = Field(SUVW, (u,v,w)->(pi)*cospi(u/2)*cospi(v/2)*cospi(w/2))

R  = 


uu =  R11*u + R12*v + R13*w
vv =  R21*u + R22*v + R23*w
ww =  R31*u + R32*v + R33*w

Dww, Dvv, Duu = derivativetransform(SUVW, uu, vv, ww)

#--------------------------------------------------------------------
# Set boundary conditions
#--------------------------------------------------------------------
ρ = 0 
s = exp(-((uu^2)/0.1) - ((vv^2)/0.1) - ((ww^2)/0.1)) 
b = 𝔹*𝕤

#--------------------------------------------------------------------
# Construct the wave operator in curved spacetime
#--------------------------------------------------------------------
guu = Field(SUV, (u,v)-> 0)
guv = Field(SUV, (u,v)->-1)
guw = Field(SUV, (u,v)->-1)
gvv = Field(SUV, (u,v)-> 0)
gvw = Field(SUV, (u,v)->-1)
gww = Field(SUV, (u,v)-> 0)

(tguu, tguv, tguw, tgtvv, tgvw, tgww) = inversemetrictransform(guu, guv, guw, gvv, gvw, gww, uu, vv, ww) 
g   = [tguu  tguv  tguw;
       tguv  tgvv  tgvw;
       tguw  tgvw  tgww]

D   = [Duu, Dvv, Dww]
L   = sum(g[i,j]*D[i]*D[j] for i in 1:3,  j in 1:3)

#--------------------------------------------------------------------
# Solve the system [also check the condition number and eigen values]
#--------------------------------------------------------------------
w = solve(L⊙B, ρ + b) 
writevtk(w,  "../output/3D-distorted")
