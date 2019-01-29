#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2019
# Transform between Galerkin and Cardinal basis
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# 1D case
#--------------------------------------------------------------------

function basistransform(u::Field{GaussLobatto{Tag, N, max, min}})::Field{Chebyshev{Tag, N, max, min}} where {Tag, N, max, min}
    α = Field(u.space, x->0)
    for m in 0:N
        α.value[m+1] = sum((2/N)*u.value[k]*cheb(m, chebx(k,N)) for k in 1:N+1)
    end
    return Field(Chebyshev{Tag, N, max, min}, α.value)
end

#--------------------------------------------------------------------
# 2D case
#--------------------------------------------------------------------

function basistransform(u::Field{T}) where T<:ProductSpace{Chebyshev{Tag2, N2, max2, min2}, 
                                                           Chebyshev{Tag1, N1, max1, min1}} where {Tag1, Tag2, N1, N2, max1, max2, min1, min2}
    A = u.value
    α = zeros(reverse(size(u.space)))
    f = zeros(size(u.space)) 

    for m in 0:N1, j in 0:N2
        yj = chebx(j+1, N2)
        α[j+1,m+1] = sum(n == 0 || n == N2 ? (1/2)*A[m+1,n+1]*cheb(n, yj) : A[m+1,n+1]*cheb(n, yj) for n in 0:N2)
    end

    for j in 0:N2, i in 0:N1
        xi = chebx(i+1, N1)
        f[i+1,j+1] = sum(m == 0 || m == N1 ? (1/2)*α[j+1,m+1]*cheb(m, xi) : α[j+1,m+1]*cheb(m, xi) for m in 0:N1) 
    end

    return Field(ProductSpace{GaussLobatto{Tag2, N2, max2, min2}, GaussLobatto{Tag1, N1, max1, min1}}, f)
end

function basistransform(u::Field{T}) where T<:ProductSpace{GaussLobatto{Tag2, N2, max2, min2}, 
                                                           GaussLobatto{Tag1, N1, max1, min1}} where {Tag1, Tag2, N1, N2, max1, max2, min1, min2}
    B = u.value
    α = zeros(reverse(size(u.space)))
    f = zeros(size(u.space)) 

    for m in 0:N1, j in 0:N2
        yj = chebx(j+1, N2)
        wj = chebw(j+1, N2)
        α[j+1,m+1] = 2*(N2+1)*sum(n == 0 || n == N2 ? (1/2)*B[m+1,n+1]*wj*cheb(n, yj) : B[m+1,n+1]*wj*cheb(n, yj) for n in 0:N2)
    end

    for j in 0:N2, i in 0:N1
        xi = chebx(i+1, N1)
        wi = chebw(i+1, N1)
        f[i+1,j+1] = 2*(N1+1)*sum(m == 0 || m == N1 ? (1/2)*α[j+1,m+1]*wi*cheb(m, xi) : α[j+1,m+1]*wi*cheb(m, xi) for m in 0:N1) 
    end
    return Field(ProductSpace{Chebyshev{Tag2, N2, max2, min2}, Chebyshev{Tag1, N1, max1, min1}}, f)
end

function interpolate(u::Field{S}, ::Type{T}) where {S<:ProductSpace{GaussLobatto{Tag2, N2, max2, min2}, 
                                                                    GaussLobatto{Tag1, N1, max1, min1}},
                                                    T<:ProductSpace{GaussLobatto{Tag2, N2new, max2, min2},
                                                                    GaussLobatto{Tag1, N1new, max1, min1}}} where {Tag1, Tag2, N1, N2, 
                                                                                                                   max1, max2, min1, min2, 
                                                                                                                   N2new, N1new}
    A = basistransform(u).value
    α = zeros(N2new+1, N1+1)
    f = zeros(N1new+1, N2new+1) 

    @show A
    @show size(α)
    @show size(f)

    for m in 0:N1, j in 0:N2new
        yj = chebx(j+1, N2new)
        @show yj
        α[j+1,m+1] = sum(n == 0 || n == N2 ? (1/2)*A[m+1,n+1]*cheb(n, yj) : A[m+1,n+1]*cheb(n, yj) for n in 0:N2)
    end

    for j in 0:N2new, i in 0:N1new
        xi = chebx(i+1, N1new)
        @show xi
        f[i+1,j+1] = sum(m == 0 || m == N1 ? (1/2)*α[j+1,m+1]*cheb(m, xi) : α[j+1,m+1]*cheb(m, xi) for m in 0:N1) 
    end

    return Field(ProductSpace{GaussLobatto{Tag2, N2new, max2, min2}, GaussLobatto{Tag1, N1new, max1, min1}}, f)
end
