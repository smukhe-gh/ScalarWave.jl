#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 03-2019
# Test 3D spaces
#--------------------------------------------------------------------

struct W end

SW = GaussLobatto{W, 11, 2, -3} 
SV = GaussLobatto{V, 11, 2, -3} 
SU = GaussLobatto{U, 11, 4, -5} 

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

@testset "3Dspace" begin
    @test D3W*q ≈ t
    @test D3V*q ≈ s
    @test D3U*q ≈ r
    @test D3W*f ≈ Field(SUVW, (U,V,W)->4*(W^3))
    @test D3V*f ≈ Field(SUVW, (U,V,W)->3*(V^2))
    @test D3U*f ≈ Field(SUVW, (U,V,W)->2*(U))
end
