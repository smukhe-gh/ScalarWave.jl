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
@test vandermonde(2*N, x)[4, 12] == cheb(11, x[4])
@test size(vandermonde(2*N,x)) == (10,2*N+1)

# test Gaussian integration
function exactint(N::Int)::Float64
    c = randn(N)
    p = rand(0:10,N)
    a = rand(-100:0,1)
    b = rand(1:100,1)

    # just a smart way to construct the polynomial and
    # it's integral
    f(x) = sum(c.*(x.^p))
    intf(x) = sum((c./(p+1)).*(x.^(p+1)))
    numericintf = quadgk(f(x), a, b; abstol=0, maxevals=10^7, order=2*N, norm=vecnorm)
    exactintf = intf(b) - intf(a) 
    @show numericintf[2]
    @show numericintf[1] - exactintf
    return numericintf[1] - exactintf 
end

@test_broken exactint(10) < 1e-14
