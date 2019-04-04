#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 03-2019
# Testing NLSolve.jl
#--------------------------------------------------------------------

using NLsolve

#--------------------------------------------------------------------
# Try a 1D problem
#--------------------------------------------------------------------

S = GaussLobatto(U, 26) 
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

