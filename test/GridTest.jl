#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

@test 1==1
using QuadGK

function hconvergence(N::Int)::Float64
    dbase = distribute(N, 1, x-> sin(pi*x), y->sin(pi*y))
    chebGrid        = Float64[chebx(i, N) for i in 1:N+1]
    chebGridData    = dbase[[1,1]]
    gaussGrid, w    = gauss(2*N)  
    exactGridData   = Float64[sin(pi*i) + sin(pi*j) for i in gaussGrid, j in gaussGrid]
    interpGridData  = interpolatePatch(chebGridData, gaussGrid, gaussGrid)   
    errorGridData   = interpGridData - exactGridData
    L2errorGridData = sqrt((w'*(errorGridData.^2)*w)/(w'*(exactGridData.^2)*w))
    return L2errorGridData
end

