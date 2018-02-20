#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

function projectboundary(func::Function, N::Int)::Array{Float64,1}
    coeffs = zeros(N+1)
    for m in 0:N
        if m == 0
            # TODO: Test guassian integration routines
            coeffs = quadgk(x->func(cos(x))*cos(m*x), 0, pi; abstol=0, maxevals=10^7, order=2*N, norm=vecnorm)/pi
        else
            coeffs = quadgk(x->func(cos(x))*cos(m*x), 0, pi; abstol=0, maxevals=10^7, order=2*N, norm=vecnorm)/(pi/2)
        end
    end
    # TODO: Make this operation work element-wise. Can you put this inside the upper loop?
    return vandermonde(N, chebgrid(N))*coeffs
end

function interpolatePatch(patch::Patch, x::Array{Float64,1}, y::Array{Float64,1})::Patch
    N      = size(patch.value)[1] - 1
    fmodal = extractPatchCoeffs(patch)
    fnodal = zeros(size(x)[1], size(y)[1])
    # TODO: Test loop. 
    for m in 1:N+1, j in 1:N+1
        elem = 0.0
        for k in 1:N+1, j in 1:N+1
            elem = elem + chebx[m,j]*fmodal[j,k]*chebx[n,k]
        end
    fnodal[m,n] = elem
    end
    #fnodal = vandermonde(N,x)*fmodal*vandermonde(N,y)'
    return Patch(patch.loc, fnodal)
end


