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

function operatorNBC{T<:Int}(N::T, M::T)::Array{Float64, 4}
    operator = zeros(N+1, N+1, N+1, N+1)
    for index in CartesianRange(size(operator))
        k = index.I[1]
        i = index.I[2]
        l = index.I[3]
        j = index.I[4]
        operator[index] = (2.0/M)*chebw(i,N)*chebw(k,N)*chebd(k,l,N)*chebd(i,j,N)
    end
    return operator
end

function operatorNBCW{T<:Int}(N::T, M::T)::Array{Float64, 4}
    operator = zeros(N+1, N+1, N+1, N+1)
    for index in CartesianRange(size(operator))
        k = index.I[1]
        i = index.I[2]
        l = index.I[3]
        j = index.I[4]
        operator[index] = (2.0/M)*chebd(k,l,N)*chebd(i,j,N)
    end
    return operator
end

function check_getIC(N::Int, M::Int)::Array{Float64,1}
    y = Float64[(chebx(i,N) - 1)/M for i in 1:N+1] 
    x = Float64[(chebx(j,N) + 1)/M for j in 1:N+1]
    return y
end
@test check_getIC(4,2) == getIC(4, 2, [2,1], :C)


