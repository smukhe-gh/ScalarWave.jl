#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------
# TODO: One can test the convergence routines separately, just by representing functions on grids.

function pconvergence(N::Int, fbnd1:: Function, fbnd2::Function, fsol::Function)::Float64
    dbase = distribute(N, 1, fbnd1, fbnd2)
    chebGrid        = Float64[chebx(i, N) for i in 1:N+1]
    chebGridData    = dbase[[1,1]]
    gaussGrid, w    = gauss(2*N)
    gaussGrid       = - gaussGrid
    exactGridData   = Float64[fsol(i,j) for i in gaussGrid, j in gaussGrid]
    interpGridData  = interpolatePatch(chebGridData, gaussGrid, gaussGrid).value
    errorGridData   = interpGridData - exactGridData
    L2errorGridData = sqrt((w'*(errorGridData.^2)*w)/(w'*(exactGridData.^2)*w))
    return L2errorGridData
end

function hconvergence(N::Int, M::Int, fbnd1:: Function, fbnd2::Function, fsol::Function)::Float64
    dbase = distribute(N, M, fbnd1, fbnd2)
    chebGrid        = Float64[chebx(i, N) for i in 1:N+1]
    gaussGrid, w    = gauss(2*N)
    gaussGrid       = -gaussGrid
    errorvec        = zeros(M*M)
    for i in 1:M, j in 1:M
        loc             = [i,j]
        gaussLocalx     = Float64[coordtransL2G(M, loc[1], x) for x in gaussGrid]
        gaussLocaly     = Float64[coordtransL2G(M, loc[1], y) for y in gaussGrid] 
        exactGridData   = Float64[fsol(x,y) for x in gaussLocalx, y in gaussLocaly]
        interpGridData  = interpolatePatch(dbase[loc], gaussGrid, gaussGrid).value
        errorGridData   = interpGridData - exactGridData
        L2errorGridData = sqrt((w'*(errorGridData.^2)*w)/(w'*(exactGridData.^2)*w))
        errorvec[i+j]   = L2errorGridData
    end
    return sum(errorvec)
end
