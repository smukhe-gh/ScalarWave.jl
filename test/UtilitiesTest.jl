#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

function testcoordtrans(N::Int, M::Int)
    xp = chebgrid(N)
    xp2xg1 = Float64[coordtransL2G(M, 1, chebx(i, N)) for i in 1:N+1]
    xp2xg2 = Float64[coordtransL2G(M, 2, chebx(i, N)) for i in 1:N+1]
    xg2xp  = Float64[coordtransG2L(M, 2, xg) for xg in xp2xg2]
    (xp2xg1[end] == xp2xg2[1]) && (xg2xp ≈ xp)
end

function testshapeH2L2H(a4N::Array{Float64,4}, Nx::Int, Ny::Int)::Bool
    shapeL2H(shapeH2L(a4N), Nx-1, Ny-1) == a4N
end

function testshapeH2L2H(b2N::Array{Float64,2}, Nx::Int, Ny::Int)::Bool
    shapeL2H(shapeH2L(b2N), Nx-1, Ny-1) == b2N
end

function testchebint(Nx::Int, Ny::Int)::Float64
    X  = chebgrid(Nx)
    Y  = chebgrid(Ny)
    B  = Float64[x^2 - y^3 + x^4*y^2 for x in X, y in Y]
    wx = chebweights(Nx)
    wy = chebweights(Ny)
    numint = wx'*B*wy
    return abs.(numint - 8/5)
end

function testL1norm(Nx::Int, Ny::Int)::Float64
    X = chebgrid(Nx)
    Y = chebgrid(Ny)
    A = Float64[x^2 - y^3 for x in X, y in Y]
    B = Float64[x^2 - y^3 + x^4*y^2 for x in X, y in Y]
    errorL1 = L1norm(B, A, chebweights(Nx), chebweights(Ny))
    exactL1 = 4/15 
    return abs(errorL1 - exactL1)
end

function testL2norm(Nx::Int, Ny::Int)::Float64
    X = chebgrid(Nx)
    Y = chebgrid(Ny)
    A = Float64[x^2 - y^3 for x in X, y in Y]
    B = Float64[x^2 - y^3 + x^4*y^2 for x in X, y in Y]
    errorL2 = L2norm(B, A, chebweights(Nx), chebweights(Ny))
    exactL2 = 2/sqrt(45)
    return abs(errorL2 - exactL2)
end

function testarray2dict2array(Nx, Ny, M)::Bool
    sPatch = rand((Nx+1)*M, (Ny+1)*M)
    return dict2array(array2dict(sPatch, Nx, Ny, M)) == sPatch
end

function testsaveandloadgrid(Nx::Int, Ny::Int, M::Int)::Bool
    sPatch = rand((Nx+1)*M, (Ny+1)*M)
    dbase  = array2dict(sPatch, Nx, Ny, M) 
    datetime = savegrid(dbase, ".")
    ldbase = loadgrid("./$datetime.jld")
    run(`rm ./$datetime.jld`)
    return ldbase["patch"][[1,1]].value == dbase[[1,1]].value
end

@test coordtransL2G(2, 1, -1.0) ≈ 0.0 
@test coordtransL2G(2, 2,  1.0)  ≈ 0.0
@test testcoordtrans(8,4) == true
@test jacobian(12) == 1/12
@test testshapeH2L2H(randn(4,12,4,12), 4, 12) == true
@test testshapeH2L2H(randn(6,10), 6, 10) == true
@test testchebint(12,14)  < 1e-14
@test testL1norm(14, 11)  < 1e-14
@test testL2norm(11, 14)  < 1e-14
@test testarray2dict2array(4,5,3) == true
@test testsaveandloadgrid(2, 2, 2) == true
