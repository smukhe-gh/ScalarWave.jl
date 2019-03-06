#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2018
# Test 1D Newton Solver. 
# TODO: Test your solver with NLSolve.jl
#--------------------------------------------------------------------

S = GaussLobatto(U, 36) x = Field(S, x->sin(x))
n = Field(S, x->(x+1)*(x-1))
abstol  = 1e-9
maxiter = 40

function F(u::Field{S})::Field{S} where {S}
    D = derivative(S)
    return D*D*u + u #- exp(u)
end

function J(u::Field{S}, Δu::Union{Symbol, Field{S}})::Operator{S} where {S}
    D = derivative(S)
    return D*D*Δu + eye(S)*Δu #- exp(u)*Δu
end

function Fvec(space::S, uvec::Array{Float64,1})::Array{Float64,1} where {S}
    return vec(F(shape(space, uvec))) 
end

function Jvec(space::S, uvec::Array{Float64,1})::Array{Float64,2} where {S}
    return vec(J(shape(space, uvec), :u))
end

function Bvec(space::S, uvec::Array{Float64,1})::Array{Float64,2} where {S}
    return vec(boundary(space))
end

xsolved = Newton(S, vec(x+n), 
                 Fvec, Jvec, Bvec, 
                 maxiter, abstol) 

