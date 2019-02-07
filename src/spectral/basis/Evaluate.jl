#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2019
# Evaluate a field function at an arbitrary point it's space
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Start with the 1D case 
#--------------------------------------------------------------------

function Base. in(coordinate::Number, S)::Bool 
    if maximum(S) >= coordinate >= minimum(S)
        return true
    else
        return false
    end
end

function (c::Field{S})(x::Number)::Number where {S<:Chebyshev{Tag, N, max, min}} where {Tag, N, max, min}
    @assert x in c.space
    xmap = (x - (max + min)/2)*(2/(min - max))
    return sum(prefactor(order+1, N)*cheb(order, xmap)*c.value[order+1] for order in 0:N)
end

function (u::Field{S})(x::Number)::Number where {S<:GaussLobatto{Tag, N}} where {Tag, N}
    @assert x in u.space
    c = basistransform(u)
    return c(x)
end

#--------------------------------------------------------------------
# Now, work with the 2D case 
#--------------------------------------------------------------------

function (α::Field{S})(UU::Number, VV::Number)::Number where S<:ProductSpace{Chebyshev{TagV, NV, maxV, minV}, 
                                                                             Chebyshev{TagU, NU, maxU, minU}} where {TagV, TagU, NV, NU, maxV, maxU, minV, minU}
                                                       
    @assert VV in Chebyshev{TagV, NV, maxV, minV}
    @assert UU in Chebyshev{TagU, NU, maxU, minU}

    Umap = (UU - (maxU + minU)/2)*(2/(minU - maxU))
    Vmap = (VV - (maxV + minV)/2)*(2/(minV - maxV))

    (M, N) = order(α.space)
    return sum(prefactor(m+1,n+1, M, N)*α.value[m+1,n+1]*cheb(m, Umap)*cheb(n, Vmap) for m in 0:M, n in 0:N) 
end

function (u::Field{S})(UU::Number, VV::Number)::Number where S<:ProductSpace{GaussLobatto{TagV, NV, maxV, minV}, 
                                                                             GaussLobatto{TagU, NU, maxU, minU}} where {TagV, TagU, NV, NU, maxV, maxU, minV, minU}
    @assert UU in Chebyshev{TagU, NU, maxU, minU}
    @assert VV in Chebyshev{TagV, NV, maxV, minV}
    c = basistransform(u)
    return c(UU, VV)
end
