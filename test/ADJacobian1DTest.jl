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
        J[index] = Δf(index.I[1], index.I[2], x)
    end
    x = Array{eltype(x)}(x) 
end

function Δf(i::Int, j::Int, x::Array{T,1})::Number where {N, T <: Union{X, Dual{X}}} where {X}
    x[j] = Dual(x[j], 1)
    δf   = reshape(F(reshape(S, x))) 
    x[j] = realpart(x[j])
    return dualpart(δf[i])
end

#--------------------------------------------------------------------
# Auxilliary functions 
#--------------------------------------------------------------------

function Base. exp(u::Field{S})::Field{S} where {S}
    return Field(u.space, exp.(u.value))
end

function F(u::Field{S})::Field{S} where {S}
    f =  D*D*u - exp(u)
    return (I-B)*f + B*(u-b)
end

function J(u::Field{S})::Operator{S} where {S}
    return D*D - exp(u)*I
end

#--------------------------------------------------------------------
# Test solver 
#--------------------------------------------------------------------

struct M end
S = ChebyshevGL{M, 22, Float64}(-1,1)
D = derivative(S)
B = incomingboundary(S) + outgoingboundary(S)
I = identity(S)

u = Field(S, x->1)
b = Field(S, x->0)  

u = nlsolve(f!, reshape(u); autodiff=:forward, show_trace=true, ftol=1e-12)

using PyPlot
plot(reshape(S, u.zero))
show()
