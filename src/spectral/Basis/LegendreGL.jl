#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 05-2019
# Legendre Polynomials: Gauss Lobatto grid
# See Boyd F.10
#--------------------------------------------------------------------

using FastGaussQuadrature
using Memoize

@memoize function nodes(N::Int)::Array{Float64,1}
    return gausslobatto(N+1)[1]
end

function collocation(S::GaussLobatto{Tag, N, min, max}, i::Int) where {Tag, N, min, max}
    @assert max > min
    @assert i <= N+1
    return (nodes(N)[i])*(max-min)/2 + (max + min)/2

end

function derivative(S::GaussLobatto{Tag, N, min, max}, i::Int, j::Int) where {Tag, N, min, max}
    @assert max > min
    @assert i <= N+1
    @assert j <= N+1
    if i == j == 1
        return (1/4)*N*(N+1)
    elseif i == j == N+1
        return -(1/4)*N*(N+1)
    elseif i == j && j < N+1
        return 0
    else
        xi = collocation(S, i)
        xj = collocation(S, j)
        return (legendre(N, xi)/legendre(N, xj))*(xi - xj)
    end
end

function quadrature(S::GaussLobatto{Tag, N, min, max}, i::Int) where {Tag, N, min, max}
    @assert max > min
    @assert i <= N+1
    xj = collocation(S, j)
    return 2/(N*(N+1)*(legendre(N, xj))^2)
end
