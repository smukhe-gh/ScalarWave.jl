#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 10-2019
# Test basis change functions
#--------------------------------------------------------------------

export basistransform

function cheb(m::Int, x::Float64)::Float64
    if abs(x) <= 1
        return cos(m*acos(x))
    elseif x >= 1
        return cosh(m*acosh(x))
    else
        return ((-1)^m)*cosh(m*acosh(-x))
    end
end

function chebx(i::T, N::T)::Float64 where {T<:Int}
    return -cospi((i-1)/N)
end

function prefactor(i::Int, N::Int)::Number
    (i == 1 || i == N+1) ? (return 1/2) : (return 1)
end

function prefactor(i, j, N, M)
    return prefactor(i, N)*prefactor(j,M)
end

#--------------------------------------------------------------------
# 1D case
#--------------------------------------------------------------------

function basistransform(α::Field{Chebyshev{Tag, M, T}})::Field{ChebyshevGL{Tag, M, T}} where {Tag, M, T}
    v = Field(ChebyshevGL{Tag, M, T}(α.space.min, α.space.max), similar(α.value))
    N = order(v.space)
    for gridindex in 1:N+1
        v.value[gridindex] = sum(prefactor(order+1, N)*cheb(order, chebx(gridindex, N))*α.value[order+1] for order in 0:N)
    end
    return v
end

function basistransform(u::Field{ChebyshevGL{Tag, M, T}})::Field{Chebyshev{Tag, M, T}} where {Tag, M, T}
    α = Field(Chebyshev{Tag, M, T}(u.space.min, u.space.max), similar(u.value))
    N = order(u.space)
    for order in 0:N
        α.value[order+1] = (2/N)*sum(prefactor(gridindex, N)*cheb(order, chebx(gridindex, N))*u.value[gridindex] for gridindex in 1:N+1)
    end
    return α
end

function Base. in(coordinate::Number, space::S)::Bool where {S<:Space}
    if space.max >= coordinate >= space.min
        return true
    else
        return false
    end
end

function (c::Field{S})(x::Number)::Number where {S<:Chebyshev{Tag, M, T}} where {Tag, M, T}
    @assert x in c.space
    (min, max) = (c.space.min, c.space.max)
    xmap = (x - (max + min)/2)*(2/(min - max))
    N = M - 1 
    return sum(prefactor(order+1, N)*cheb(order, xmap)*c.value[order+1] for order in 0:N)
end

function (u::Field{S})(x::Number)::Number where {S<:ChebyshevGL{Tag, M, T}} where {Tag, M, T}
    @assert x in u.space
    c = basistransform(u)
    return c(x)
end

# #--------------------------------------------------------------------
# # 2D case
# #--------------------------------------------------------------------

# function basistransform(α::Field{T}) where T<:ProductSpace{Chebyshev{Tag1, M1, T}, 
                                                           # Chebyshev{Tag2, M2, T}} where {Tag1, Tag2, M1, M2, T}
    # (M, N) = (M1 - 1, M2 - 1)
    # f = zeros(M1, M2)

    # for i in 1:M+1, j in 1:N+1
        # f[i,j] = sum(prefactor(m+1,n+1, M, N)*α.value[m+1,n+1]*cheb(m, chebx(i, M))*cheb(n, chebx(j,N)) for m in 0:M, n in 0:N) 
    # end

    # return Field(ProductSpace{ChebyshevGL{Tag1, M1, T}(α.space.S1.min, α.space.S1.max), 
                              # ChebyshevGL{Tag2, M2, T}(α.space.S2.min, α.space.S2.max)}, f)
# end

# function basistransform(u::Field{T}) where T<:ProductSpace{GaussLobatto{Tag1, M1, T}, 
                                                           # GaussLobatto{Tag2, M2, T}} where {Tag1, Tag2, M1, M2, T}
    # (M, N) = order(u.space)
    # c = zeros(size(u.space))

    # for m in 0:M, n in 0:N
        # c[m+1,n+1] = (4/(M*N))*sum(prefactor(i,j, M, N)*u.value[i,j]*cheb(m, chebx(i, M))*cheb(n, chebx(j,N)) for i in 1:M+1, j in 1:N+1) 
    # end

    # return Field(ProductSpace{Chebyshev{Tag1, M1, T}(u.space.S1.min, u.space.S1.max), 
                              # Chebyshev{Tag2, M2, T}(u.space.S2.min, u.space.S2.max)}, c)
# end

# function (α::Field{S})(UU::Number, VV::Number)::Number where S<:ProductSpace{Chebyshev{TagV, NV, maxV, minV}, 
                                                                             # Chebyshev{TagU, NU, maxU, minU}} where {TagV, TagU, NV, NU, maxV, maxU, minV, minU}
                                                       
    # @assert VV in Chebyshev{TagV, NV, maxV, minV}
    # @assert UU in Chebyshev{TagU, NU, maxU, minU}

    # Umap = (UU - (maxU + minU)/2)*(2/(minU - maxU))
    # Vmap = (VV - (maxV + minV)/2)*(2/(minV - maxV))

    # (M, N) = order(α.space)
    # return sum(prefactor(m+1,n+1, M, N)*α.value[m+1,n+1]*cheb(m, Umap)*cheb(n, Vmap) for m in 0:M, n in 0:N) 
# end

# function (u::Field{S})(UU::Number, VV::Number)::Number where S<:ProductSpace{GaussLobatto{TagV, NV, maxV, minV}, 
                                                                             # GaussLobatto{TagU, NU, maxU, minU}} where {TagV, TagU, NV, NU, maxV, maxU, minV, minU}
    # @assert UU in Chebyshev{TagU, NU, maxU, minU}
    # @assert VV in Chebyshev{TagV, NV, maxV, minV}
    # c = basistransform(u)
    # return c(UU, VV)
# end
