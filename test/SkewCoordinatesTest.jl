#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 03-2019
# Test operations with skew coordinates
# Test Laplace equation with interesting boundary conditions
#--------------------------------------------------------------------

struct X end
struct Y end
struct Z end

SXY = ProductSpace{GaussLobatto{Y, 20,  1.0, -1.0},
                   GaussLobatto{X, 20,  1.0, -1.0}}
DY, DX = derivative(SXY)
B = boundary(Spacelike, SXY)
L = DX*DX + DY*DY
u = Field(SXY, (X,Y)->X*Y*(X+1)*(Y+1))
u = solve(L+B, B*u)
XF = Field(SXY, (X,Y)->X)
YF = Field(SXY, (X,Y)->Y)

α  = 0 
β  = π/24 
UF = Field(SXY, (X,Y)->cos(α)*X + cos(β)*Y)
VF = Field(SXY, (X,Y)->sin(α)*X + sin(β)*Y)

g = [    1     cos(α - β);
     cos(α - β)     1    ]
D = [DY, DX]

L = sum(g[i,j]*D[i]*D[j] for i in 1:2, j in 1:2)
u = Field(SXY, (X,Y)->X*Y*(X+1)*(Y+1), VF, UF)
u = solve(L+B, B*u)

# Now test Laplace equation in 3D
SXYZ = ProductSpace{GaussLobatto{Z, 10,  1.0, -1.0},
                    GaussLobatto{Y, 10,  1.0, -1.0},
                    GaussLobatto{X, 10,  1.0, -1.0}}

DZ, DY, DX = derivative(SXYZ)
g = [1 0 0;
     0 1 0;
     0 0 1]
D = [DZ, DY, DX]
B = boundary(Null, SXYZ)
@time L = sum(g[i,j]*D[i]*D[j] for i in 1:3, j in 1:3)
u = Field(SXYZ, (X,Y,Z)->X+Y+Z) 
u = solve(L+B, B*u)

CZ = Field(GaussLobatto{Z, 10,  1.0, -1.0}, X->X)
CY = Field(GaussLobatto{Y, 10,  1.0, -1.0}, X->X)
CX = Field(GaussLobatto{X, 10,  1.0, -1.0}, X->X)

using WriteVTK
vtkfile = vtk_grid("testVTK", CX.value, CY.value, CZ.value) 
vtk_point_data(vtkfile, u.value, "scalar")
outfiles = vtk_save(vtkfile)
