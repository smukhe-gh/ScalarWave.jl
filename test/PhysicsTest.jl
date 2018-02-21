#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

function check_operator(N::Int, M::Int)::Array{Float64,2}
    w = Float64[(1/M^2)*chebw(i,N)*chebw(j,N) for i in 1:N+1, j in 1:N+1]
    d = Float64[chebd(i,j,N) for i in 1:N+1, j in 1:N+1] 
    i = eye(N+1,N+1)
    
    W  = diagm(vec(w))
    DU = (1/M)*kron(i, d) 
    DV = (1/M)*kron(d, i)
    WD = W*(DU*DV + DV*DU)
    IM = eye((N+1)^2, (N+1)^2) 

    b    = zeros(N+1, N+1)
    b[1,:] = b[:,1] = 1 
    bvec = vec(b)
    for j in 1:1+N, i in 1:N+1 
        if b[i,j] == 1.0
            K = i + (N+1)*(j-1)    
            WD[K,:] = IM[K,:]
        end
    end
    return WD
end

@test check_operator(2, 1) ≈ reshapeA(operator(2, 1))
@test check_operator(4, 3) ≈ reshapeA(operator(4, 3))

@test getIC(4, 2, [1,1], (x)->x , :R).value == (chebgrid(4) + 1)/2
@test getIC(4, 2, [1,1], (x)->x , :C).value == (chebgrid(4) + 1)/2
@test getIC(4, 2, [1,2], (x)->x , :R).value == (chebgrid(4) + 1)/2
@test getIC(4, 2, [1,2], (x)->x , :C).value == (chebgrid(4) - 1)/2
@test getIC(4, 2, [2,1], (x)->x , :R).value == (chebgrid(4) - 1)/2
@test getIC(4, 2, [2,1], (x)->x , :C).value == (chebgrid(4) + 1)/2
@test getIC(4, 2, [2,2], (x)->x , :R).value == (chebgrid(4) - 1)/2
@test getIC(4, 2, [2,2], (x)->x , :C).value == (chebgrid(4) - 1)/2
