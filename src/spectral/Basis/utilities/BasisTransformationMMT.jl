#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2019
# Transform between Galerkin and Cardinal basis
#--------------------------------------------------------------------

function prefactor(i, N)
    (i == 1 || i == N+1) ? (return 1/2) : (return 1)
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
#--------------------------------------------------------------------

function basistransform(u::Field{T}) where T<:ProductSpace{Chebyshev{Tag2, N2, max2, min2}, 
                                                           Chebyshev{Tag1, N1, max1, min1}} where {Tag1, Tag2, N1, N2, max1, max2, min1, min2}
    A = u.value
    α = zeros(reverse(size(u.space)))
    f = zeros(size(u.space)) 

    for m in 1:N1+1, j in 1:N2+1
        α[j, m] = sum(prefactor(n+1, N2)*A[m, n+1]*cheb(n, chebx(j, N2)) for n in 0:N2)
    end

    for j in 1:N2+1, i in 1:N1+1
        f[i, j] = sum(prefactor(m+1, N1)*α[j, m+1]*cheb(m, chebx(i, N1)) for m in 0:N1) 
    end

    return Field(ProductSpace{GaussLobatto{Tag2, N2, max2, min2}, GaussLobatto{Tag1, N1, max1, min1}}, f)
end

function basistransform(u::Field{T}) where T<:ProductSpace{GaussLobatto{Tag2, N2, max2, min2}, 
                                                           GaussLobatto{Tag1, N1, max1, min1}} where {Tag1, Tag2, N1, N2, max1, max2, min1, min2}
    A = u.value
    α = zeros(reverse(size(u.space)))
    β = zeros(size(u.space)) 

    for i in 0:N1, n in 0:N2
        α[n+1, i+1] = sum(prefactor(j, N2)*A[i+1, j]*cheb(n, chebx(j, N2)) for j in 1:N2+1)
    end

    for n in 0:N2, m in 0:N1
        β[m, n] = sum(prefactor(i, N1)*α[n+1, i]*cheb(m, chebx(i, N1)) for i in 1:N1+1) 
    end

    return Field(ProductSpace{Chebyshev{Tag2, N2, max2, min2}, Chebyshev{Tag1, N1, max1, min1}}, β)
end
