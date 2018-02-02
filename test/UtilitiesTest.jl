#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

@test coordtrans(2, [0.0,0.0], [1,2]) ≈ [-0.5,0.5]
@test coordtrans(2, [0.0,0.0], [2,1]) ≈ [0.5,-0.5]

function check_coordtrans(N::Int)::Bool
    x11 = zeros(N+1)
    y11 = zeros(N+1) 
    x21 = zeros(N+1)
    y21 = zeros(N+1)
    x12 = zeros(N+1)
    y12 = zeros(N+1)
    
    for i in 1:N+1
        x11[i], y11[i] = coordtrans(2, [chebx(i,N),chebx(i,N)], [1,1])
        x21[i], y21[i] = coordtrans(2, [chebx(i,N),chebx(i,N)], [2,1])
        x12[i], y12[i] = coordtrans(2, [chebx(i,N),chebx(i,N)], [1,2])
    end

    X = vcat(x11, x12)
    Y = vcat(y11, y21)
    if X == Y
        return true
    else
        return false
    end
end
@test check_coordtrans(4) == true    

function check_shape_reshapeA(a4N::Array{Float64,4})::Bool
    shapeA(reshapeA(a4N)) == a4N ? true : false
end
@test check_shape_reshapeA(randn(4,4,4,4)) == true

function check_shape_reshapeB(b2N::Array{Float64,2})::Bool
    shapeB(reshapeB(b2N)) == b2N ? true : false
end 
@test check_shape_reshapeB(randn(6,6)) == true

N = 9
x = Float64[chebx(i, N) for i in 1:N+1]
@test vandermonde(N, x)[3,4] == cheb(2, x[4])
@test size(vandermonde(N,x)) == (10,10)
@test size(vandermonde(19,x)) == (10,20)
