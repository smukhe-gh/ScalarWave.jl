#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 03-2019
# Test 3D spaces
#--------------------------------------------------------------------

struct W end

SW = GaussLobatto{W, 14, 2, -3} 
SV = GaussLobatto{V, 12, 3, -4} 
SU = GaussLobatto{U, 13, 4, -5} 

SUV  = ProductSpace{SV, SU}
SUVW = ProductSpace{SW, SV, SU}


# test operations with fields
f = Field(SUVW, (U,V,W)->U^2 + V^3 + W^4)
g = Field(SUVW, (U,V,W)->U*V*W)
h = Field(SUVW, (U,V,W)->U^2 + V^3 + W^4 + U*V*W)
u = Field(SUVW, (U,V,W)->(U*V*W)*(U^2 + V^3 + W^4))

@testset "Group operations" begin
    @test f ≈ f
    @test f + g ≈ h
    @test f*g ≈ u
end

l = Field(SUV, (U,V)->U*V)
m = Field(SUV, (U,V)->V)
n = Field(SUV, (U,V)->U)
DV, DU = derivative(SUV)

@testset "2Dspace" begin
    @test DV*l ≈ n 
    @test DU*l ≈ m 
    @test vec(DV) == kron(derivative(SV).value, eye(SU).value) 
    @test vec(DU) == kron(eye(SV).value, derivative(SU).value)
end

D3W, D3V, D3U = derivative(SUVW)

q = Field(SUVW, (U,V,W)->U*V*W)
t = Field(SUVW, (U,V,W)->U*V)
s = Field(SUVW, (U,V,W)->U*W)
r = Field(SUVW, (U,V,W)->V*W)
w = Field(SUVW, (U,V,W)->W)
z = Field(SUVW, rand(size(SUVW)...))

@testset "3Dspace" begin
    @test reshape(SUVW, reshape(z)) == z
    @test reshape(SUVW, reshape(D3U)).value == D3U.value
    @test reshape(SUVW, reshape(D3U*D3V*D3W)).value ==  (D3U*D3V*D3W).value
    @test w*D3W*q ≈ q
    @test D3U*D3V*q ≈ w
    @test D3V*q ≈ s
    @test D3U*q ≈ r
    @test D3W*f ≈ Field(SUVW, (U,V,W)->4*(W^3))
    @test D3V*f ≈ Field(SUVW, (U,V,W)->3*(V^2))
    @test D3U*f ≈ Field(SUVW, (U,V,W)->2*(U))
end

#--------------------------------------------------------------------
# Now test Laplace equation in 3D
#--------------------------------------------------------------------

using LinearAlgebra

struct Z end
struct Y end
struct X end

SZ = GaussLobatto{Z, 12, 1, -5}
SY = GaussLobatto{Y, 11, 1, -5}
SX = GaussLobatto{X, 12, 1, -5}
SXYZ = ProductSpace{SZ, SY, SX}


DZ, DY, DX = derivative(SXYZ)
B = boundary(Spacelike, SXYZ)
L = DZ*DZ + DY*DY + DX*DX
b = Field(SXYZ, (X,Y,Z)->X+Y+Z) 
u = solve(L⊙B, B*b)

@show maximum(abs(B*u - B*b))
@show maximum(abs(u))

@show maximum(abs(u.value[1, :, :] - u.value[end, :, :]))
@show maximum(abs(u.value[:, 1, :] - u.value[:, end, :]))
@show maximum(abs(u.value[:, :, 1] - u.value[:, :, end]))

writevtk(u, "../output/laplace")
writevtk(B*u, "../output/laplace-boundary")


