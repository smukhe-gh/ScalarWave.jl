#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

function interpolatePatch(patch::Patch, N::Int)::Patch
    ON      = size(patch.value)[1] - 1
    vndm   = vandermonde(ON, chebgrid(N))
    fmodal = extractPatchCoeffs(patch)
    """
    fnodal = zeros(N+1, N+1)
    for i in 1:N+1, j in 1:N+1
        elem = 0.0
        for m in 1:N+1, n in 1:N+1
            elem = elem + vndm[i,m]*fmodal[m,n]*vndm[j,n]
        end 
        fnodal[i,j] = elem
    end
    """
    fnodal = vndm*fmodal*vndm' 
    return Patch(patch.loc, fnodal)
end

function projectboundary(func::Function, N::Int)::Array{Float64,1}
    coeffs = zeros(N+1)
    for m in 0:N
        if m == 0
            coeffs = quadgk(x->func(cos(x))*cos(m*x), 0, pi; abstol=0, maxevals=10^7, order=2*N, norm=vecnorm)/pi
        else
            coeffs = quadgk(x->func(cos(x))*cos(m*x), 0, pi; abstol=0, maxevals=10^7, order=2*N, norm=vecnorm)/(pi/2)
        end
    end
    return vandermonde(N, chebgrid(N))*coeffs
end



