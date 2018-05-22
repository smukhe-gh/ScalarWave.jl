function computexpansion(grid::Grid, metric::Metric, index::CartesianIndex)::Float64
    # For spherically symmetric spacetimes, one can simplify the
    # numerical calculation of the expansion. See Eq. (7.24)
    # in Baumgarte & Shapiro (2010)
    r = grid.r_of_UV[index]
    M = grid.params.mass 
    alpha = (1 - 2M/r)/(1 + 2M/r)
end


function schwarzschildOP(grid::Grid, metric::Metric)::Array{Float64, 4}
    @assert grid.params.size[1] == grid.params.size[2]
    N = grid.params.size[1]
	operator = zeros(N, N, N, N)
    for index in CartesianRange(size(operator))    
        k = index.I[1]
        i = index.I[2]
        l = index.I[3]
        j = index.I[4]   
        operator[index] = 2*chebw(i,Ny)*chebw(k,Nx)*chebd(k,l,Nx)*chebd(i,j,Ny)
                            - dOmegaU[k,l]/Omega[k,l] chebd(i,k,Nx) 
                            - dOmegaV[k,l]/Omega[k,l] chebd(i,k,Nx)
	end
	return operator
end

function RHS{T<:Int}(fn::Function, Nx::T, Ny::T, M::T, loc::Array{Int,1})::Array{Float64,2}
    @assert Nx == Ny
    rhs = projectonPatchbyRestriction(fn, Nx, Ny, M, loc)
    for index in CartesianRange(size(rhs))
        i = index.I[1]
        j = index.I[2]
        rhs[i,j] = chebw(i,Nx)*chebw(j,Ny)*rhs[i,j]
    end
    return rhs
end
