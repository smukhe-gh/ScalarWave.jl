#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

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

@test testderivOP(2, 2) ≈ shapeH2L(derivOP(2, 2))
@test_broken testderivOP(2, 4) ≈ shapeH2L(derivOP(2, 4))
@test testboundaryOP(2,2) == shapeH2L(boundaryOP(2,2))
@test testboundaryOP(2,4) == shapeH2L(boundaryOP(2,4))
@test RHS((x,y)->x^8+y^9, 8, 4, 7, [2,1]) ≈ testRHS(8, 4, 7, [2,1])
