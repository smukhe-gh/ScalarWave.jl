#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2019
# Transform between Galerkin and Cardinal basis
#--------------------------------------------------------------------

function prefactor(i, N)
    (i == 1 || i == N+1) ? (return 1/2) : (return 1)
end

function prefactor(i, j, N, M)
    return prefactor(i, N)*prefactor(j,M)
end


#--------------------------------------------------------------------
# 1D case
#--------------------------------------------------------------------

function basistransform(α::Field{Chebyshev{Tag, N, max, min}})::Field{GaussLobatto{Tag, N, max, min}} where {Tag, N, max, min}
    u = Field(GaussLobatto{Tag, N, max, min})
    for gridindex in 1:N+1
        u.value[gridindex] = sum(prefactor(order+1, N)*cheb(order, chebx(gridindex, N))*α.value[order+1] for order in 0:N)
    end
    return u
end

function basistransform(u::Field{GaussLobatto{Tag, N, max, min}})::Field{Chebyshev{Tag, N, max, min}} where {Tag, N, max, min}
    α = Field(Chebyshev{Tag, N, max, min})
    for order in 0:N
        α.value[order+1] = (2/N)*sum(prefactor(gridindex, N)*cheb(order, chebx(gridindex, N))*u.value[gridindex] for gridindex in 1:N+1)
    end
    return α
end

#--------------------------------------------------------------------
# 2D case
# TODO: Use MMT to make the transforms faster.
#--------------------------------------------------------------------

function basistransform(α::Field{T}) where T<:ProductSpace{Chebyshev{Tag1, N1, max1, min1}, 
                                                           Chebyshev{Tag2, N2, max2, min2}} where {Tag1, Tag2, N1, N2, max1, max2, min1, min2}
    (M, N) = order(α.space)
    f = zeros(size(α.space))

    for i in 1:M+1, j in 1:N+1
        f[i,j] = sum(prefactor(m+1,n+1, M, N)*α.value[m+1,n+1]*cheb(m, chebx(i, M))*cheb(n, chebx(j,N)) for m in 0:M, n in 0:N) 
    end
    return Field(ProductSpace{GaussLobatto{Tag1, N1, max1, min1}, 
                              GaussLobatto{Tag2, N2, max2, min2}}, f)
end

function basistransform(u::Field{T}) where T<:ProductSpace{GaussLobatto{Tag1, N1, max1, min1}, 
                                                           GaussLobatto{Tag2, N2, max2, min2}} where {Tag1, Tag2, N1, N2, max1, max2, min1, min2}
    (M, N) = order(u.space)
    c = zeros(size(u.space))

    for m in 0:M, n in 0:N
        c[m+1,n+1] = (4/(M*N))*sum(prefactor(i,j, M, N)*u.value[i,j]*cheb(m, chebx(i, M))*cheb(n, chebx(j,N)) for i in 1:M+1, j in 1:N+1) 
    end

    return Field(ProductSpace{Chebyshev{Tag1, N1, max1, min1}, 
                              Chebyshev{Tag2, N2, max2, min2}}, c) 
end