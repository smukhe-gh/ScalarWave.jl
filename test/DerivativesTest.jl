#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 05-2018
#--------------------------------------------------------------------

function testindexdifferentiation(N::Int)
    fx = Float64[x^2 + y^3 + x^3*y^2 for x in chebgrid(N), y in chebgrid(N)]
    D  = Float64[chebd(i,j,N) for i in 1:N+1, j in 1:N+1]
    ndfdx = zeros(N+1, N+1)
    ndfdy = zeros(N+1, N+1)
    nddfdx = zeros(N+1, N+1)
    nddfdy = zeros(N+1, N+1)
    nddfdxdy = zeros(N+1, N+1)

    for index in CartesianRange((N+1, N+1))
        i = index[1]
        j = index[2]

        ndfdy[i,j]    = sum(chebd(j,n,N)*fx[i,n] for n in 1:N+1)
        ndfdx[i,j]    = sum(chebd(i,m,N)*fx[m,j] for m in 1:N+1)
        nddfdy[i,j]   = sum(chebd(j,k,N)*chebd(k,n,N)*fx[i,n] for n in 1:N+1, k in 1:N+1)
        nddfdx[i,j]   = sum(chebd(i,l,N)*chebd(l,m,N)*fx[m,j] for m in 1:N+1, l in 1:N+1)
        nddfdxdy[i,j] = sum(chebd(j,n,N)*chebd(i,k,N)*fx[k,n] for k in 1:N+1, n in 1:N+1)
    end

    @test ndfdy ≈ Float64[3*y^2 + 2*y*x^3 for x in chebgrid(N), y in chebgrid(N)]
    @test ndfdx ≈ Float64[2*x + 3*y^2*x^2 for x in chebgrid(N), y in chebgrid(N)]
    @test nddfdy ≈ Float64[6*y + 2*x^3    for x in chebgrid(N), y in chebgrid(N)]
    @test nddfdx ≈ Float64[2 + 6*y^2*x    for x in chebgrid(N), y in chebgrid(N)]
    @test nddfdxdy ≈ Float64[6*y*x^2  for x in chebgrid(N), y in chebgrid(N)]
end

function testoperator(Nx::Int, Ny::Int)::Array{Float64,4}
    """
    The main operator we're using 
    sum(chebw(i,Nx)*chebw(j,Ny)*
         chebd(i,m,Nx)*delta(j,n)*phi[m,n]*
         delta(i,p)*chebd(j,q,Ny)*phi[p,q] 
         for m in 1:Nx+1, n in 1:Ny+1, p in 1:Nx+1, q in 1:Ny+1, 
         i in 1:Nx+1, j in 1:Ny+1)
    where we sum over i, j in the computation of the operator.
    
    And for computing the action of the operator on the vector field
    operator[m,n,p,q] phi[m,n] phi[p,q]
    and when you reshape it into a 2D array (m,n) and (p,q) are bunched 
    together.
    """

    operator = zeros(Nx+1, Ny+1, Nx+1, Ny+1)
    for index in CartesianRange(size(operator)) 
        (m,n,p,q) = index.I
        operator[m,n,p,q] = chebw(p,Nx)*chebw(n,Ny)*chebd(p,m,Nx)*chebd(n,q,Ny)
    end
    return operator
end

#---------------------------------------------------
# Testing region
#---------------------------------------------------

Nx, Ny = (14,32)
phi = Float64[sin(x-y)*exp(x+y) 
              for x in chebgrid(Nx), y in chebgrid(Ny)] 
operator = testoperator(Nx, Ny)
integral = vec(phi)'*shapeH2L(operator)*vec(phi)
@test integral ≈ -6.990469114220025
