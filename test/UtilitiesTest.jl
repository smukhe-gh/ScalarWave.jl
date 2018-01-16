#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

@test coordtrans(2, [0.0,0.0], [1,2]) ≈ [0.5,-0.5]
@test coordtrans(2, [0.0,0.0], [2,1]) ≈ [-0.5,0.5]

function check_coordtrans(N::Int, M::Int)
    # FIXME: This function isn't implemented properly. 
    #        However, the coordtrans function works.
    px = zeros(M, N+1)
    for m in 1:M, i in 1:N+1
        px[m,i] = coordtrans(M, [chebx(i,N), 1.0], [m,1])[1]
        if m != M && i==N+1
            px[m,i] == -100.0 
        end
    end
    px = reshape(px, (N+1)*M)
    fpx = px[px .!= -100.0]
    return fpx    
end
@test_broken check_coordtrans(3,2) ≈ Float64[chebx(i, 7) for i in 1:7] 

function check_shape_reshapeA(a4N::Array{Float64,4})::Bool
    shapeA(reshapeA(a4N)) == a4N ? true : false
end
@test check_shape_reshapeA(randn(4,4,4,4)) == true

function check_shape_reshapeB(b2N::Array{Float64,2})::Bool
    shapeB(reshapeB(b2N)) == b2N ? true : false
end 
@test check_shape_reshapeB(randn(6,6)) == true

function check_computeB(N,M)
# construct Pascal's triangle
end
