#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 07-2019
# Modal to nodal basis transformations
#--------------------------------------------------------------------

function cheb(m::Int, x::T)::T where {T}
    if abs(x) <= 1
        return cos(m*acos(x))
    elseif x >= 1
        return cosh(m*acosh(x))
    else
        return ((-1)^m)*cosh(m*acosh(-x))
    end
end

function inverseM(space::S)::AbstractArray{T,2} where {S<:Chebyshev{Tag, N, T}} where {Tag, N, T}
    M = zeros(T, (N, N))
    for index in CartesianIndices(M)
        gridindex, orderplus1 = index.I
        M[index] = cheb(orderplus1 - 1, collocation(ChebyshevGL{Tag, N, T}(space.min, space.max), gridindex))
    end
end

function nodal2modal(u::Field{ChebyshevGL{Tag, N, T}})::Field{Chebyshev{Tag, N, T}} where {Tag, N, T}
    return Field(Chebyshev{Tag, N, T}(u.space.min, u.space.max), inv(inverseM(N,T))*u.value) 
end

function modal2nodal(a::Field{Chebyshev{Tag, N, T}})::Field{ChebyshevGL{Tag, N, T}} where {Tag, N, T}
    return Field(ChebyshevGL{Tag, N, T}(a.space.min, a.space.max), inverseM(N,T)*a.value) 
end

struct M end
S = ChebyshevGL{M, 4, Float64}(-1, 1)
f = Field(S, x->x^3)
@show nodal2modal(f) 
