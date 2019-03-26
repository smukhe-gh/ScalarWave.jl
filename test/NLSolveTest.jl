#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 03-2019
# Testing NLSolve.jl
#--------------------------------------------------------------------

using NLsolve

#--------------------------------------------------------------------
# Try a 1D problem
#--------------------------------------------------------------------

S = GaussLobatto(U, 46) 
x = Field(S, x->0)

function F(u::Field{S})::Field{S} where {S}
    D = derivative(S)
    return D*D*u - exp(u)
end

function f(x::Array{Float64,1})::Array{Float64,1}
    return vec(F(Field(S, x)))
end

w = nlsolve(f, x.value)
@show minimum(w.zero)

#=
#--------------------------------------------------------------------
# Try a 2D linear problem
#--------------------------------------------------------------------

SUV = ProductSpace{GaussLobatto(V, 20),
                   GaussLobatto(U, 20)}
x = Field(SUV, (U,V)->U*V)

function F(u::Field{SUV})::Field{SUV} where {SUV}
    DV, DU = derivative(SUV)
    return DU*DU*u + DV*DV*u 
end

function f(x::Array{Float64,1})::Array{Float64,1}
    return vec(F(Field(SUV, shape(SUV, x))))
end

# w = nlsolve(f, vec(x.value))
# contourf(Field(SUV, shape(SUV, w.zero)), 100, globallevels=20)
# show()

#--------------------------------------------------------------------
# Try a 2D linear problem
#--------------------------------------------------------------------
# XXX: Add the linear boundary operator 

SUV = ProductSpace{GaussLobatto(V, 20),
                   GaussLobatto(U, 20)}
x = Field(SUV, (U,V)->U*V)
using LinearAlgebra

function F(u::Field{SUV})::Field{SUV} where {SUV}
    DV, DU = derivative(SUV)
    B = boundary(Spacelike, SUV)
    @show minimum(eigvals(vec((DU*DU + DV*DV + B))))
    return DU*DU*u + DV*DV*u + (B*u  - x) 
end

function f(x::Array{Float64,1})::Array{Float64,1}
    return vec(F(Field(SUV, shape(SUV, x))))
end

w = nlsolve(f, vec(x.value))
contourf(Field(SUV, shape(SUV, w.zero)), 100, globallevels=20)
show()

# TODO: Replace rows to impose boundary conditions
function ⊕(A::ProductSpaceOperator{S}, B::ProductSpaceOperator{S})::ProductSpaceOperator{S} where {S}
    D = vec(A)
    E = vec(B)
    for rowindex in prod(size(S)):
        if sum(E[rowindex]) > 0
            D[rowindex] = E[rowindex]
        end
    end
    return ProductSpaceOperator(S, reshape(D, (size(S)..., size(S)...)))
end
=#
