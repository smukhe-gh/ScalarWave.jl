#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2018
# Test 2D Newton Solver. 
#--------------------------------------------------------------------

# function F(ϕ::Field{S})::Field{S} where {S}
    # DV, DU = derivative(S)
    # return DU*DV*ϕ
# end

# function J(u::Field{S}, Δϕ::Union{Symbol, Field{S}})::ProductSpaceOperator{S} where {S}
    # DV, DU = derivative(S)
    # return (DU*ϕ)*(DV*Δϕ) + (DV*ϕ)*(DU*Δϕ)
# end

function F(ϕ::Field{S})::Field{S} where {S}
    DV, DU = derivative(S)
    return DU*DU*ϕ + DV*DV*ϕ
end

function J(u::Field{S}, Δϕ::Union{Symbol, Field{S}})::ProductSpaceOperator{S} where {S}
    DV, DU = derivative(S)
    return (DV*DV*Δϕ) + (DU*DU*Δϕ)
end

function B(space::Type{S}, sym::Symbol)::ProductSpaceOperator{S} where {S} 
    bnd = boundary(Spacelike, space)
    (S == :Δ0) ? (return bnd - bnd) : (return bnd)
end

function Fvec(space::Type{S}, Svec::Array{Float64,1})::Array{Float64,1} where {S}
    ϕ = Sshape(space, Svec)
    return vec(F(ϕ)) 
end

function Jvec(space::Type{S}, Svec::Array{Float64,1})::Array{Float64,2} where {S}
    ϕ = Sshape(space, Svec)
    return vec(J(ϕ, :Δϕ))             
end

function Bvec(space::Type{S}, Svec::Array{Float64,1})::Array{Float64,2} where {S}
    return vec(B(space, :Δϕ))
end

function Sshape(space::Type{S}, Svec::Array{Float64,1}) where {S}
    return Field(space, shape(space, Svec)) 
end

function Svec(ϕ::Field{S})::Array{Float64,1} where {S}
    return vec(ϕ)
end

#--------------------------------------------------------------------
# Solve scalar wave on Minkowski
# FIXME: The solver doesn't converge in one iteration
#--------------------------------------------------------------------

SUV = ProductSpace{GaussLobatto(V, 21),
                   GaussLobatto(U, 21)}
SUVW = ProductSpace{SUV, GaussLobatto(W, 21)}
abstol  = 1e-9
maxiter = 40

noise = Field(SUV, (U,V)->(V-3)*(V-1)*(U+1)*(U-1))
ϕ = Field(SUV, (U,V)->U+V)

xsolved = Newton(SUV, vec(ϕ + noise), 
                 Fvec, Jvec, Bvec, 
                 maxiter, abstol) 


#--------------------------------------------------------------------
# Solve Schwarzschild 
#--------------------------------------------------------------------

SUV = ProductSpace{GaussLobatto{V, 21,  3,  1},
                   GaussLobatto{U, 21, -3, -5}}

M = 1.0
abstol = 1e-7
maxiter = 40

# Initial Data for Minkowski 
r = Field(SUV, (U,V)->((V-U)/2))
f = Field(SUV, (U,V)->1)
ϕ = Field(SUV, (U,V)->0)

# Initial Data for Schwarzschild + noise
n = Field(SUV, (U,V)->10^-4*(U+3)*(U+5)*(V-3)*(V-1))
r = Field(SUV, (U,V)->find_r_of_UV(U,V,M))
f = (16*M^3/r)*exp(-r/2M)
ϕ = Field(SUV, (U,V)->0)

xsolved = Newton(SUV, Svec(f, r+n, ϕ), 
                 Fvec, Jvec, Bvec, 
                 maxiter, abstol) 

