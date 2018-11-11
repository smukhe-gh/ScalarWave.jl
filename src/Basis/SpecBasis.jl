#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
# Wrapper functions for basis choices
#--------------------------------------------------------------------

function collocation(space::Type{S}, i::Int)::Rational where  {S<:Taylor{Tag, N, max, min}} where {Tag, N, max, min}
    @assert max > min
    @assert i <= N+1
    @assert typeof(min - max) <: Rational
    return -(-1 + (2*(i-1)//N))*(min - max)/2 + (max + min)/2
end

function collocation(space::Type{S}, i::Int)::Float64 where {S<:GaussLobatto{Tag, N, max, min}} where {Tag, N, max, min}
    @assert max > min
    @assert i <= N+1
    return chebx(i, N)*(min - max)/2 + (max + min)/2
end

function derivative(space::Type{S}, i::Int, j::Int)::Rational where {S<:Taylor{Tag, N, max, min}} where {Tag, N, max, min}
    @assert i <= N+1
    @assert j <= N+1
    x = [collocation(space, k) for k in 1:N+1]
    if i == j
        return sum((k==j ? 0 : 1/(x[j] - x[k])) for k in 1:N+1)*(2/(min - max))
    else
        ai = prod((k==i ? 1 : (x[i] - x[k])) for k in 1:N+1)
        aj = prod((k==j ? 1 : (x[j] - x[k])) for k in 1:N+1)
        return ai/(aj*(x[i] - x[j]))*(2/(min - max))
    end
end

function derivative(space::Type{S}, i::Int, j::Int)::Float64 where {S<:GaussLobatto{Tag, N, max, min}} where {Tag, N, max, min}
    @assert i <= N+1
    @assert j <= N+1
    return chebd(i, j, N)*(2/(min - max))
end

function derivative(space::Type{S}, i::Int)::Float64 where {S<:GaussLobatto{Tag, N, max, min}} where {Tag, N, max, min}
    @assert i <= N+1
    return chebw(i, N)*(2/(min - max))
end
