#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 06-2018
#--------------------------------------------------------------------

function RHS{T<:Int}(fn::Function, Nx::T, Ny::T, M::T, loc::Array{Int,1})::Array{Float64,2}
    rhs = projectonPatchbyRestriction(fn, Nx, Ny, M, loc)
    for index in CartesianRange(size(rhs))
        (i,j) = index.I
        rhs[i,j] = chebw(i,Nx)*chebw(j,Ny)*rhs[i,j]
    end
    return rhs
end
