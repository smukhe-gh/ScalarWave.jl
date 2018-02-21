#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

function check_operator(Nx::Int, Ny::Int)::Array{Float64,2}
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
   
    WD = zeros((Nx+1)*(Ny+1), (Nx+1)*(Ny+1))
    IM = eye((Nx+1)*(Ny+1), (Nx+1)*(Ny+1)) 
    b  = zeros(Nx+1, Ny+1)
    b[1,:] = b[:,1] = 1 
    for ny in 1:Ny+1, nx in 1:Nx+1 
        if b[nx,ny] == 1.0
            K = sub2ind((Nx+1, Ny+1), nx, ny)
            WD[K,:] = IM[K,:]
        end
    end
    return WD
end

@test check_operator(2, 2) ≈ reshapeA(operator(2, 2))
@test check_operator(4, 2) ≈ reshapeA(operator(4, 2))

"""
@test getIC(4, 2, [1,1], (x)->x , :R).value == (chebgrid(4) + 1)/2
@test getIC(4, 2, [1,1], (x)->x , :C).value == (chebgrid(4) + 1)/2
@test getIC(4, 2, [1,2], (x)->x , :R).value == (chebgrid(4) + 1)/2
@test getIC(4, 2, [1,2], (x)->x , :C).value == (chebgrid(4) - 1)/2
@test getIC(4, 2, [2,1], (x)->x , :R).value == (chebgrid(4) - 1)/2
@test getIC(4, 2, [2,1], (x)->x , :C).value == (chebgrid(4) + 1)/2
@test getIC(4, 2, [2,2], (x)->x , :R).value == (chebgrid(4) - 1)/2
@test getIC(4, 2, [2,2], (x)->x , :C).value == (chebgrid(4) - 1)/2
"""
