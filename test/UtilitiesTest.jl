#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

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

function check_shape_reshapeA(a4N::Array{Float64,4})::Bool
    shapeA(reshapeA(a4N)) == a4N ? true : false
end

function check_shape_reshapeB(b2N::Array{Float64,2})::Bool
    shapeB(reshapeB(b2N)) == b2N ? true : false
end 

function check_quadgk(N::Int)::Float64
    a, b = (-1,1)
    f(x) = 16*x^2 - 6*x + 2
    intf(x) = (16/3)*x^3 - 3*x^2 + 2*x
    exactint = intf(b) - intf(a)
    numint = quadgk(f,a,b; reltol=1e-15, order=4)
    return abs.(exactint - numint[1])
end

function check_chebint(N::Int)::Float64
    X = Y = chebgrid(N)    
    B = Float64[x^2 - y^3 + x^4*y^2 for x in X, y in Y]
    w = chebweights(N)
    numint = w'*B*w
    return abs.(numint - 8/5)
end

function check_gaussint(N::Int)::Float64
    x, w = gauss(N) 
    X = Y = -x
    B = Float64[x^2 - y^3 + x^4*y^2 for x in X, y in Y]
    numint = w'*B*w
    return abs.(numint - 8/5)
end

function check_L1norm(N::Int)::Float64
    X = Y = chebgrid(N)
    A = Float64[x^2 - y^3 for x in X, y in Y]
    B = Float64[x^2 - y^3 + x^4*y^2 for x in X, y in Y]
    errorL1 = L1norm(B, A, chebweights(N))
    exactL1 = 1/5 
    return abs(errorL1 - exactL1)
end

function check_L2norm(N::Int)::Float64
    X = Y = chebgrid(N)
    A = Float64[x^2 - y^3 for x in X, y in Y]
    B = Float64[x^2 - y^3 + x^4*y^2 for x in X, y in Y]
    errorL2 = L2norm(B, A, chebweights(N))
    exactL2 = sqrt(7/3)/6
    return abs(errorL2 - exactL2)
end

function gaussint2D(N::Int)::Float64
    x,w = gauss(N)
    xp = (x-1)/2
    yp = (x+1)/2
    fx = Float64[x^2 - y^3 + x^4*y^2 for x in xp, y in yp]
    J = jacobian(2)
    numint = J*w'*fx*w
    return abs(numint - 3/20) 
end

function check_savegrid()
    dbase = distribute(4, 2, x->sin(pi*x), y->sin(pi*y))
    savegrid(dbase, "../output")
    return true
end

@test check_savegrid() == true
@test check_coordtrans(4) == true    
@test jacobian(12) == (1/12)^2
@test check_shape_reshapeA(randn(4,4,4,4)) == true
@test check_shape_reshapeB(randn(6,6)) == true
@test coordtrans(2, [0.0,0.0], [1,2]) ≈ [-0.5,0.5]
@test coordtrans(2, [0.0,0.0], [2,1]) ≈ [0.5,-0.5]
@test vandermonde(2*9, chebgrid(9))[4, 12] == cheb(11, chebgrid(9)[4])
@test size(vandermonde(2*9,chebgrid(9))) == (10,2*9+1)
@test check_chebint(10) < 1e-14
@test check_gaussint(8) < 1e-14
@test check_quadgk(12) < 1e-14
@test check_L1norm(14) < 1e-14
@test check_L2norm(11) < 1e-14
@test gaussint2D(10) < 1e-14
