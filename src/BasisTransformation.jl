#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 09-2018
# Basis transformation functions with Type-I DFT
#--------------------------------------------------------------------

#---------------------------------------------------------------
# 1D functions
#---------------------------------------------------------------

function basistransform(u::Field{GaussLobatto{Tag, N}}, ::Type{Chebyshev{Tag, N}})::Field{Chebyshev{Tag, N}} where {Tag, N}  
    c = (1/N)*(FFTW.r2r(u.value, FFTW.REDFT00))
    return Field(Chebyshev{Tag, N}, c)
end

function basistransform(u::Field{Chebyshev{Tag, N}}, ::Type{GaussLobatto{Tag, N}})::Field{GaussLobatto{Tag, N}} where {Tag, N}  
    n = (FFTW.r2r(N*u.value, FFTW.REDFT00))/(2*N)
    return Field(GaussLobatto{Tag, N}, n)
end

function interpolate(u::Field{GaussLobatto{Tag, N1}}, x::Field{GaussLobatto{Tag, N2}})::Field{GaussLobatto{Tag, N2}} where {Tag, N1, N2}
    c = basistransform(u, Chebyshev{Tag, N1})
    S = GaussLobatto{Tag, N2}
    v = zeros(N2+1)
    for index in 1:N2+1
        v[index] = sum(m==0 || m==N1 ? (c.value[m+1]/2)*cheb(m, x.value[index]) : c.value[m+1]*cheb(m, x.value[index]) for m in 0:N1) 
    end
    return Field(GaussLobatto{Tag, N2}, v)
end

#---------------------------------------------------------------
# 2D functions
#---------------------------------------------------------------

function basistransform(u::Field{ProductSpace{GaussLobatto{Tag, N1}, 
                                              GaussLobatto{Tag, N2}}}, 
                          ::Type{ProductSpace{Chebyshev{Tag, N1},
                                              Chebyshev{Tag, N2}}})::Field{ProductSpace{Chebyshev{Tag, N1},
                                                                                        Chebyshev{Tag, N2}}} where {Tag, N1, N2}  
    c = (1/(N1*N2))*(FFTW.r2r(u.value, FFTW.REDFT00))
    return Field(ProductSpace{Chebyshev{Tag, N1}, Chebyshev{Tag, N2}}, c)
end

function basistransform(u::Field{ProductSpace{Chebyshev{Tag, N1}, 
                                              Chebyshev{Tag, N2}}}, 
                          ::Type{ProductSpace{GaussLobatto{Tag, N1},
                                              GaussLobatto{Tag, N2}}})::Field{ProductSpace{GaussLobatto{Tag, N1},
                                                                                           GaussLobatto{Tag, N2}}} where {Tag, N1, N2}  
    n = (FFTW.r2r((N1*N2)*u.value, FFTW.REDFT00))/(2*(N1*N2))
    return Field(ProductSpace{GaussLobatto{Tag, N1}, GaussLobatto{Tag, N2}}, n)
end

function interpolate(u::Field{ProductSpace{GaussLobatto{Tag, N1}, 
                                           GaussLobatto{Tag, N2}}}, 
                       ::Type{ProductSpace{GaussLobatto{Tag, N3},
                                           GaussLobatto{Tag, N4}}})::Field{ProductSpace{Chebyshev{Tag, N3},
                                                                                        Chebyshev{Tag, N4}}} where {Tag, N1, N2, N3, N4}  
    c  = basistransform(u, ProductSpace{Chebyshev{Tag, N1}, Chebyshev{Tag, N2}})
    @show size(c)
    PS = ProductSpace{GaussLobatto{Tag, N2}, GaussLobatto{Tag, N2}}
    
    α  = zeros(size(PS))
    for m in 0:N4, j in 0:N3
        α[m+1, j+1] = sum( (n==0 || n==N3) ? (c[m+1,n+1]/2)*cheb(n,chebx(j,N3)) : c[m+1,n+1]*cheb(n,chebx(j,N3)) for n in 0:N3)
    end

    𝑓  = zeros(size(PS))
    for m in 0:N4, j in 0:N3
        𝑓[m+1, j+1] = sum( (m==0 || m==N4) ? (α[m+1,j+1]/2)*cheb(m,chebx(j,N4)) : α[m+1,j+1]*cheb(m,chebx(j,N4)) for m in 0:N4)
    end
    
    return Field(PS, 𝑓)
end

struct M end

function interpolate(u::Field, N1::Int, N2::Int, N3::Int, N4::Int) 
    c  = basistransform(u, ProductSpace{Chebyshev{M, N1}, Chebyshev{M, N2}})
    @show size(c)
    PS = ProductSpace{GaussLobatto{M, N2}, GaussLobatto{M, N2}}
    
    α  = zeros(size(PS))
    for m in 0:N4, j in 0:N3
        α[m+1, j+1] = sum( (n==0 || n==N3) ? (c[m+1,n+1]/2)*cheb(n,chebx(j,N3)) : c[m+1,n+1]*cheb(n,chebx(j,N3)) for n in 0:N3)
    end

    𝑓  = zeros(size(PS))
    for m in 0:N4, j in 0:N3
        𝑓[m+1, j+1] = sum( (m==0 || m==N4) ? (α[m+1,j+1]/2)*cheb(m,chebx(j,N4)) : α[m+1,j+1]*cheb(m,chebx(j,N4)) for m in 0:N4)
    end
    
    return Field(PS, 𝑓)
end
