#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 09-2018
# 1D and 2D Basis transformation functions with Type-I DFT
# TODO: Add ability to do this with partial summations
#--------------------------------------------------------------------

function basistransform(u::Field{GaussLobatto{Tag, N}}, ::Type{Chebyshev{Tag, N}})::Field{Chebyshev{Tag, N}} where {Tag, N}  
    c = (1/N)*(FFTW.r2r(u.value, FFTW.REDFT00))
    return Field(Chebyshev{Tag, N}, c)
end

function basistransform(u::Field{Chebyshev{Tag, N}}, ::Type{GaussLobatto{Tag, N}})::Field{GaussLobatto{Tag, N}} where {Tag, N}  
    n = (FFTW.r2r(N*u.value, FFTW.REDFT00))/(2*N)
    return Field(GaussLobatto{Tag, N}, n)
end

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
