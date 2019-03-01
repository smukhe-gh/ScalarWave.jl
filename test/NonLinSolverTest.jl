#--------------------------------------------------------------------
# Test non-linear solver
# Soham M 03/2019
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# The non-linear and linearized equations
#--------------------------------------------------------------------

# S = GaussLobatto(U, 36)
# x = Field(S, x->0)
# abstol  = 1e-10
# maxiter = 40

# function F(u::Field{S})::Field{S} where {S}
    # D = derivative(S)
    # return D*D*u - exp(u)
# end

# function J(u::Field{S}, Δu::Union{Symbol, Field{S}})::Operator{S} where {S}
    # D = derivative(S)
    # return D*D*Δu - exp(u)*Δu
# end

# function Fvec(space::S, uvec::Array{Float64,1})::Array{Float64,1} where {S}
    # return vec(F(shape(space, uvec))) 
# end

# function Jvec(space::S, uvec::Array{Float64,1})::Array{Float64,2} where {S}
    # return vec(J(shape(space, uvec), :u))
# end

# function Bvec(space::S, uvec::Array{Float64,1})::Array{Float64,2} where {S}
    # return vec(boundary(space))
# end

# xsolved = Newton(S, vec(x), 
                 # Fvec, Jvec, Bvec, 
                 # maxiter, abstol) 

#--------------------------------------------------------------------
# Interface with the Newton solver 
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

