#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

function testindexdifferentiation(N::Int)
    fx = Float64[x^2 + y^3 + x^3*y^2 for x in chebgrid(N), y in chebgrid(N)]

    D  = Float64[chebd(i,j,N) for i in 1:N+1, j in 1:N+1]
    dfdy = kron(D, eye(N+1))*vec(fx)
    dfdx = kron(eye(N+1), D)*vec(fx)
    ddfdy = kron(D*D, eye(N+1))*vec(fx)
    ddfdx = kron(eye(N+1), D*D)*vec(fx)
    ddfdxdy = kron(eye(N+1), D)*kron(D,eye(N+1))*vec(fx)
    ndfdx = zeros(N+1, N+1)
    ndfdy = zeros(N+1, N+1)
    nddfdx = zeros(N+1, N+1)
    nddfdy = zeros(N+1, N+1)
    nddfdxdy = zeros(N+1, N+1)

    for index in CartesianRange((N+1, N+1))
        i = index[1]
        j = index[2]

        ndfdy[i,j]    = sum(delta(i,m)*chebd(j,n,N)*fx[m,n] for m in 1:N+1, n in 1:N+1)
        ndfdx[i,j]    = sum(chebd(i,m,N)*delta(j,n)*fx[m,n] for m in 1:N+1, n in 1:N+1)
        nddfdy[i,j]   = sum(delta(i,m)*chebd(j,k,N)*chebd(k,n,N)*fx[m,n] for m in 1:N+1, n in 1:N+1, k in 1:N+1)
        nddfdx[i,j]   = sum(chebd(i,l,N)*chebd(l,m,N)*delta(j,n)*fx[m,n] for m in 1:N+1, n in 1:N+1, l in 1:N+1)
        nddfdxdy[i,j] = sum(delta(i,m)*chebd(j,n,N)*chebd(m,k,N)*delta(n,l)*fx[k,l] for m in 1:N+1, n in 1:N+1, k in 1:N+1, l in 1:N+1)

        ndfdy[i,j]    = sum(chebd(j,n,N)*fx[i,n] for n in 1:N+1)
        ndfdx[i,j]    = sum(chebd(i,m,N)*fx[m,j] for m in 1:N+1)
        nddfdy[i,j]   = sum(chebd(j,k,N)*chebd(k,n,N)*fx[i,n] for n in 1:N+1, k in 1:N+1)
        nddfdx[i,j]   = sum(chebd(i,l,N)*chebd(l,m,N)*fx[m,j] for m in 1:N+1, l in 1:N+1)
        nddfdxdy[i,j] = sum(chebd(j,n,N)*chebd(i,k,N)*fx[k,n] for k in 1:N+1, n in 1:N+1)
    end
   
    @test dfdx ≈ vec(Float64[2*x + 3*y^2*x^2 for x in chebgrid(N), y in chebgrid(N)])
    @test ddfdx ≈ vec(Float64[2 + 6*y^2*x    for x in chebgrid(N), y in chebgrid(N)])
    @test dfdy ≈ vec(Float64[3*y^2 + 2*y*x^3 for x in chebgrid(N), y in chebgrid(N)])
    @test ddfdy ≈ vec(Float64[6*y + 2*x^3    for x in chebgrid(N), y in chebgrid(N)])
    @test ddfdxdy ≈ vec(Float64[6*y*x^2      for x in chebgrid(N), y in chebgrid(N)])

    @test ndfdy ≈ Float64[3*y^2 + 2*y*x^3 for x in chebgrid(N), y in chebgrid(N)]
    @test ndfdx ≈ Float64[2*x + 3*y^2*x^2 for x in chebgrid(N), y in chebgrid(N)]
    @test nddfdy ≈ Float64[6*y + 2*x^3    for x in chebgrid(N), y in chebgrid(N)]
    @test nddfdx ≈ Float64[2 + 6*y^2*x    for x in chebgrid(N), y in chebgrid(N)]
    @test nddfdxdy ≈ Float64[6*y*x^2  for x in chebgrid(N), y in chebgrid(N)]
end

function testderivOP(Nx::Int, Ny::Int)::Array{Float64,2}
    wx = Float64[chebw(i,Nx) for i in 1:Nx+1]                   #size: (Nx+1)
    wy = Float64[chebw(i,Ny) for i in 1:Ny+1]                   #size: (Ny+1)
    dx = Float64[chebd(i,j,Nx) for i in 1:Nx+1, j in 1:Nx+1]    #size: (Nx+1, Nx+1)
    dy = Float64[chebd(i,j,Ny) for i in 1:Ny+1, j in 1:Ny+1]    #size: (Ny+1, Ny+1)
    
    ix = eye(Nx+1,Nx+1)     #size: (Nx+1, Nx+1)                                     
    iy = eye(Ny+1,Ny+1)     #size: (Ny+1, Ny+1)

    w  = kron(wx,wy)        #size: (Nx+1, Ny+1)
    W  = diagm(vec(w))      #size: ((Nx+1)*(Ny+1), (Nx+1)*(Ny+1)) 
    DU = kron(iy, dx)       #size: ((Ny+1)*(Nx+1), (Ny+1)*(Nx+1)) 
    DV = kron(dy, ix)       #size: ((Ny+1)*(Nx+1), (Ny+1)*(Nx+1)) 
    WD = W*(DU*DV + DV*DU)  #size: ((Ny+1)*(Nx+1), (Ny+1)*(Nx+1)) 
    return WD
end

function testderivOPonfunc(Nx::Int, Ny::Int)
    sfunc = Float64[(sinpi(x)^2)*(cospi(x)^4) for x in chebgrid(Nx), y in chebgrid(Ny)]
    operator = testderivOP(Nx, Ny)
end

function testboundaryOP(Nx::Int, Ny::Int)::Array{Float64,2}
    bnd = zeros(Nx+1, Ny+1)
    bnd[1,:] = bnd[:,1] = 1
    return diagm(vec(bnd))
end    

function testRHS(Nx::Int, Ny::Int, M::Int, loc::Array{Int,1})
    patch = projectonPatchbyRestriction((x,y)->x^8+y^9, Nx, Ny, M, loc)
    for index in CartesianRange(size(patch))
        i = index[1]
        j = index[2]
        patch[index] = chebw(i,Nx)*chebw(j,Ny)*patch[index]
    end
    return patch
end

testindexdifferentiation(10)
@test testderivOP(2, 2) ≈ shapeH2L(derivOP(2, 2))
@test_broken testderivOP(2, 4) ≈ shapeH2L(derivOP(2, 4))
@test testboundaryOP(2,2) == shapeH2L(boundaryOP(2,2))
@test testboundaryOP(2,4) == shapeH2L(boundaryOP(2,4))
@test RHS((x,y)->x^8+y^9, 8, 4, 7, [2,1]) ≈ testRHS(8, 4, 7, [2,1])
