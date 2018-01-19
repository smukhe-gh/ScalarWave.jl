#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

function check_operator(N::Int, M::Int)::Array{Float64,2}
    w = Float64[chebw(i,N)*chebw(j,N) for i in 1:N+1, j in 1:N+1]
    d = Float64[chebd(i,j,N) for i in 1:N+1, j in 1:N+1] 
    i = eye(N+1,N+1)
    
    W  = diagm(vec(w))
    DU = kron(i, d) 
    DV = kron(d, i)
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

function check_operatorNBC(N::Int, M::Int)::Bool   
    opBNC  = operatorNBC(N,M)
    f(x,y) = x^2*y^3 + y^2*x^3
    fvec   = Float64[f(chebx(i,N), chebx(j,N)) for i in 1:N:1, j in 1:N+1]
    I      = sum(reshapeA(opBNC)*reshapeB(fvec))
    return I
end
@test_broken check_operatorNBC(10,1) ≈ 0.0

function check_initializeRHS(N::Int, M::Int)::Array{Float64,2}
    b    = zeros(N+1, N+1)
    #loc = [2,1]
    brow = Float64[(chebx(i,N) - 1)/M for i in 1:N+1] 
    bcol = Float64[(chebx(j,N) + 1)/M for j in 1:N+1]
    b[:, 1] = bcol
    b[1, :] = brow
    return b
end
@test check_initializeRHS(4,2) == initializeRHS(4, 2, [2,1], (x,y)->y, (x,y)->x)


