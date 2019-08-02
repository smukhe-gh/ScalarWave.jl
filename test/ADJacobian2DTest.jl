#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 04-2019
# Wrappers for automatic differentiation
#--------------------------------------------------------------------

using DualNumbers
using NLsolve

#--------------------------------------------------------------------
# Core functions 
#--------------------------------------------------------------------

function f!(f::Array{T,1}, x::Array{T,1}) where {T}
    f[:] = reshape(F(reshape(S, x)))
end

function j!(J::Array{T,2}, x::Array{T,1}) where {T}
    x = Array{Union{eltype(x), Dual{eltype(x)}}}(x) 
    for index in CartesianIndices(J)
        J[:, index.I[2]] = Δf(index.I[2], x)
    end
    x = Array{eltype(x)}(x) 
end

function Δf(j::Int, x::Array{T,1})::Array{T,1} where {N, T <: Union{X, Dual{X}}} where {X}
    x[j] = Dual(x[j], 1)
    δf   = reshape(F(reshape(S, x))) 
    x[j] = realpart(x[j])
    return dualpart.(δf)
end

#--------------------------------------------------------------------
# Auxilliary functions 
#--------------------------------------------------------------------

function Base. sqrt(u::Field{S})::Field{S} where {S}
    return Field(u.space, sqrt.(u.value))
end

function F(u::Field{S})::Field{S} where {S}
    return DU*DU*u + DV*DV*u + sqrt((DU*u)^2 + (DV*u)^2)
end

function J(u::Field{S})::Operator{S} where {S}
    return DU*DU + DV*DV + (1/sqrt((DU*u)^2 + (DV*u)^2))*((DU*u)*DU + (DV*u)*DV)
end

#--------------------------------------------------------------------
# Test solver 
#--------------------------------------------------------------------

struct M end
S1 = ChebyshevGL{M, 3, Float64}(1,3)
S2 = ChebyshevGL{M, 3, Float64}(1,3)
S  = ProductSpace(S1, S2)
DU, DV = derivative(S)

u = Field(S, (x,y)->x^2+y^2)
AD_J = zeros(9,9)
j!(AD_J, reshape(u))
@show maximum(abs.(reshape(J(u)) - AD_J))
display(J(u))
display(AD_J)
