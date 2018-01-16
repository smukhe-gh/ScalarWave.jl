#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

function check_operatorBC(N::Int, M::Int)::Bool
    op4N = operator(N, 1)
    op2N = reshape(op4N, ((N+1)^2,(N+1)^2))
    
    b       = zeros((N+1), (N+1))
    b[1, :] = 1.0
    b[:, 1] = 1.0
    bvec    = vec(b)
    
    for index in CartesianRange(size(b))
        K = index.I[1]
        if b[K] == 1
            if dot(op2N[K, :], op2N[K, :]) == 1
                return true
            else
                return false
            end
        end
    end    
end
@test check_operatorBC(3, 1) == true

function check_operator_regression(N::Int, M::Int)::Array{Float64,2}
    w = Float64[chebw(i,N)*chebw(j,N) for i in 1:N+1, j in 1:N+1]
    d = Float64[chebd(i,j,N) for i in 1:N+1, j in 1:N+1] 
    b = Float64[i==1 ? 1 : (j==1 ? 1 : 0) for i in 1:N+1, j in 1:N+1]
    i = eye(N+1,N+1)
    
    W  = diagm(vec(w))
    DU = kron(i, d) 
    DV = kron(d, i)
    WD = W*(DU*DV + DV*DU)
    IM = eye((N+1)^2, (N+1)^2) 

    bvec = vec(b)
    for index in CartesianRange(size(b))
        K = index.I[1]
        if b[K] == 1
            WD[K,:] = IM[K,:]
        end
    end
    return WD
end
@test_broken check_operator_regression(2, 1) ≈ reshape(operator(2, 1), (9, 9))

function check_setB_regression(N::Int)::Array{Float64,2}
    b    = zeros(N+1, N+1)
    # using M = 2, loc = [2 1]
    brow = Float64[(chebx(i,N) + 1)/2 for i in 1:N+1] 
    bcol = Float64[(chebx(j,N) - 1)/2 for j in 1:N+1]
    b[1, :] = brow
    b[:, 1] = bcol
    return b
end
@test check_setB_regression(4) == setB(4, 2, [2,1])


