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

function boundaryOP{T<:Int}(Px::T, Py::T)::Array{Float64, 4}
    bnd = zeros(Px+1, Py+1, Px+1, Py+1)
    for index in CartesianRange(size(bnd))
        i, ii, j, jj = index.I
        if  i==1 || ii==1
            bnd[index] = delta(ii,jj)*delta(i,j)
        end
    end
    return bnd
end
