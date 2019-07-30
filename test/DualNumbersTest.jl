#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 07-2019
# Experiment with Dual numbers and computing the Jacobian
#--------------------------------------------------------------------

using DualNumbers

function Base. exp(u::Field{S})::Field{S} where {S}
    return Field(u.space, exp.(u.value))
end

function J(u::Field{S})::Operator{S} where {S}
    return D*D - exp(u)*I
end

function F(u::Field{S})::Field{S} where {S}
    return D*D*u - exp(u)
end

function Δu(symbol::Symbol, i::Int, u::Field{S}; h=1e-5)::Field{S}  where {S}
    @assert i <= size(u.space)
    δu = Field(u.space, copy(u.value))
    if symbol == :+
        δu.value[i] = u.value[i] + h
    else
        δu.value[i] = u.value[i] - h
    end
    return δu
end

function ΔF(i::Int, j::Int, u::Field{S})::Number where {S}
    δuplus  = Δu(:+, j, u)
    δuminus = Δu(:-, j, u) 
    δf      = F(δuplus) - F(δuminus) 
    δu      = δuplus - δuminus
    return δf.value[i]/δu.value[j]
end

function FDJ(u::Field{S})::Operator{S} where {S}
    J = zeros(size(u.space), size(u.space))
    for index in CartesianIndices(J)
        i,j = index.I
        J[index] = ΔF(i, j, u)
    end
    return Operator(u.space, J)
end

function DualΔu(i::Int, u::Field{S}; h=1e-5)::Field{S}  where {S}
    @assert i <= size(u.space)
    δu = Field(u.space, zeros(Dual{eltype(u.value)}, size(u.space)))
    for index in CartesianIndices(δu.value)
        index.I[1] == i ?  δu.value[index] = Dual(u.value[index], 1) : δu.value[index] = Dual(u.value[index], 0)
    end
    return δu
end

function DualΔF(i::Int, j::Int, u::Field{S})::Number where {S}
    δu = DualΔu(j, u) 
    δf = F(δu)
    return δf.value[i]
end

function DualJ(u::Field{S})::Operator{S} where {S}
    J = zeros(eltype(u.value), size(u.space), size(u.space))
    for index in CartesianIndices(J)
        i,j = index.I
        J[index] = dualpart(DualΔF(i, j, u))
    end
    return Operator(u.space, J)
end


#----------------------------------------------------------
# test performance
#----------------------------------------------------------

struct M end
S = ChebyshevGL{M, 10, Float64}(-1,1)
D = derivative(S)
I = identity(S)
u = Field(S, x->sin(x)*exp(x))
@time J(u)
@time DualJ(u)


