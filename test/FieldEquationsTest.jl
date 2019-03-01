#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2019
# Test a non-linear solver with Einstein's equations 
#--------------------------------------------------------------------

SUV = ProductSpace{GaussLobatto{V, 21,  3,  1},
                   GaussLobatto{U, 21, -3, -5}}

M = 1.0
abstol = 1e-7

# Initial Data for Minkowski 
r = Field(SUV, (U,V)->((V-U)/2))
f = Field(SUV, (U,V)->1)
ϕ = Field(SUV, (U,V)->0)

# Initial Data for Schwarzschild
r = Field(SUV, (U,V)->find_r_of_UV(U,V,M))
f = (16*M^3/r)*exp(-r/2M)
ϕ = Field(SUV, (U,V)->0)

@testset "NonLin Residuals" begin
    @test norm(F(f, r, ϕ, :UU)) < abstol
    @test norm(F(f, r, ϕ, :VV)) < abstol
    @test norm(F(f, r, ϕ, :UV)) < abstol
    @test norm(F(f, r, ϕ, :θθ)) < abstol
    @test norm(F(f, r, ϕ, :TT)) < abstol
end

# Test if the linear equations and operators are working properly
f  = Field(SUV, (U,V)->(U^2)*(V^3))
r  = Field(SUV, (U,V)->(V^2)*(U^3))
ϕ  = Field(SUV, (U,V)->0)
Δr = Field(SUV, (U,V)->V^2+U^3)
Δf = Field(SUV, (U,V)->U^2+V^3)
Δϕ = Field(SUV, (U,V)->0)

JVV = Field(SUV, (U,V)->(-4/(V^5))*(3 + 2*V))  
JUU = Field(SUV, (U,V)->(-12/U^4))  
JUV = Field(SUV, (U,V)->(2/((U^7)*(V^4)))*(U^3 
                                           - (2*U^3)*(1+6*(U^3))*(V) + (V^3)*(-2 + U - 12*(U^3))))
Jθθ = Field(SUV, (U,V)->6*(-1 + U - ((U^3/V^3)*(-1+V))))

DV, DU  = derivative(SUV)
II      = eye(SUV)
ΔfOPJVV =  (2/f)*(1/r)*(DV*r)*DV - (2/f^2)*(1/r)*(DV*f)*(DV*r)*II  
ΔrOPJVV = -(2/r)*DV*DV + (2/f)*(1/r)*(DV*f)*DV - (2/f)*(1/r^2)*(DV*f)*(DV*r)*II + (2/r^2)*(DV*DV*r)*II  

@testset "Lin Residuals" begin
    @test norm(J(f, r, ϕ, :VV, Δf, Δr, Δϕ) - JVV) < abstol 
    @test norm(J(f, r, ϕ, :UU, Δf, Δr, Δϕ) - JUU) < abstol
    @test norm(J(f, r, ϕ, :θθ, Δf, Δr, Δϕ) - Jθθ) < abstol
    @test norm(J(f, r, ϕ, :UV, Δf, Δr, Δϕ) - JUV) < abstol
    @test maximum(abs.(J(f, r, ϕ, :VV, :Δf, 0, 0).value - ΔfOPJVV.value)) == 0.0
    @test maximum(abs.(J(f, r, ϕ, :VV, 0, :Δr, 0).value - ΔrOPJVV.value)) == 0.0
end

