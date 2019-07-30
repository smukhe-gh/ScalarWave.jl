#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 04-2019
# Wrappers for automatic differentiation
#--------------------------------------------------------------------

using DualNumbers

function DualNumbers. Dual(u::Field{S, D, T})::Field{S, D, Dual{T}} where {S, D, T}
    udual = Field(u.space, similar(u.value, Dual{eltype(u.value)}))
    for index in CartesianIndices(udual.value)
        udual.value[index] = Dual(u.value[index], 0)
    end
    return udual
end

function Base. exp(u::Field{S})::Field{S} where {S}
    return Field(u.space, exp.(u.value))
end

function J(u::Field{S})::Operator{S} where {S}
    return D*D - exp(u)*I
end

function F(u::Field{S})::Field{S} where {S}
    return D*D*u - exp(u)
end

function f!(F::Array{T,1}, x::Array{T,1})::Array{T,1} where {T <: Union{S, Dual{S}}} where {S}
    F[:] = reshape(resF(reshape(S, x)))
    return F
end

function Δf(index::NTuple{N, Int}, x::Array{T,1})::Number where {N, T}
    i, j = index
    x = Array{Union{Dual{eltype{x}}, eltype{x}}, x}
    x[j] = Dual(x[j], 1)
    δf   = f!(copyF, x)[i] 
    x[j] = realpart(x[j])
    return δf
end

function j!(J::Array{T,2}, x::Array{T,1})::Array{T,2} where {T}
    for index in CartesianIndices(J)
        J[index] = dualpart(Δf(index.I, x))
    end
    return J
end

struct M end
S = ChebyshevGL{M, 3, Float64}(-1,1)
D = derivative(S)
I = identity(S)
u = Field(S, x->x)
display(J(u))
display(j!(u))


# struct M end
# S1 = ChebyshevGL{M, 2, Float64}(-1, 1)
# S2 = ChebyshevGL{M, 4, Float64}(-1, 1)
# S  = ProductSpace(S1, S2)
# x  = Field(S, (x,y)->x+y)

# display(Dual(x))
# display(Dual(x, (1,2)))
