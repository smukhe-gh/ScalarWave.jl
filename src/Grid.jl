#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

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

function interpolatePatch(patch::Patch, x::Array{Float64,1}, y::Array{Float64,1})::Patch
    N      = size(patch.value)[1] - 1
    vndmx  = vandermonde(N,x)
    vndmy  = vandermonde(N,y)
    fmodal = extractPatchCoeffs(patch)
    fnodal = zeros(size(x)[1], size(y)[1])
    for i in 1:N+1, j in 1:N+1
        elem = 0.0
        for m in 1:N+1, n in 1:N+1
            elem = elem + vndmx[i,m]*fmodal[m,n]*vndmy[j,n]
        end 
        fnodal[i,j] = elem
    end
    return Patch(patch.loc, fnodal)
end


