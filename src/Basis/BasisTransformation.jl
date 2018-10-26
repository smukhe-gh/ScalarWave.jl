#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 09-2018
# 1D and 2D Basis transformation functions with Type-I DFT
# and partial summations. Note that we use a convention where the
# first and the last coefficents are divided by 2. 
# TODO: In MMT, one could use deepcopy and
#       remove the if statements [check what's faster]
#--------------------------------------------------------------------

function basistransform(u::Field{GaussLobatto{Tag, N}}, ::Type{Chebyshev{Tag, N}})::Field{Chebyshev{Tag, N}} where {Tag, N}  
    c = (1/N)*(FFTW.r2r(u.value, FFTW.REDFT00))
    return Field(Chebyshev{Tag, N}, c)
end

function basistransform(u::Field{Chebyshev{Tag, N}}, ::Type{GaussLobatto{Tag, N}})::Field{GaussLobatto{Tag, N}} where {Tag, N}  
    n = (FFTW.r2r(N*u.value, FFTW.REDFT00))/(2*N)
    return Field(GaussLobatto{Tag, N}, n)
end

function basistransform(u::Field{T}, method::Symbol) where T<:ProductSpace{GaussLobatto{Tag1, N1}, GaussLobatto{Tag2, N2}} where {Tag1, Tag2, N1, N2}
    @assert method == :dft
    c = (1/(N1*N2))*(FFTW.r2r(u.value, FFTW.REDFT00))
    return Field(ProductSpace{Chebyshev{Tag1, N1}, Chebyshev{Tag2, N2}}, c)
end

function basistransform(u::Field{T}, method::Symbol) where T<:ProductSpace{Chebyshev{Tag1, N1}, Chebyshev{Tag2, N2}} where {Tag1, Tag2, N1, N2}
    @assert method == :dft
    n = (FFTW.r2r((N1*N2)*u.value, FFTW.REDFT00))/(4*(N1*N2))
    return Field(ProductSpace{GaussLobatto{Tag1, N1}, GaussLobatto{Tag2, N2}}, n)
end

function basistransform(u::Field{T}) where T<:ProductSpace{Chebyshev{Tag2, N2}, Chebyshev{Tag1, N1}} where {Tag1, Tag2, N1, N2}
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

    return Field(ProductSpace{GaussLobatto{Tag2, N2}, GaussLobatto{Tag1, N1}}, f)
end

function basistransform(u::Field{T}) where T<:ProductSpace{GaussLobatto{Tag2, N2}, GaussLobatto{Tag1, N1}} where {Tag1, Tag2, N1, N2}
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
    return Field(ProductSpace{Chebyshev{Tag2, N2}, Chebyshev{Tag1, N1}}, f)
end
