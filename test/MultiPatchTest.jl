#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 10-2019
# Distribute computation into patches
#--------------------------------------------------------------------
using PyPlot

function contourfgrid(grid::Array{Field}, nlevels::Int)
    gridmin = minimum(grid)
    gridmax = maximum(grid)
    gridlevels = levels(grid, nlevels)
    for index in CartesianIndices(grid) 
        contourflocal(grid[index], gridlevels, gridmin, gridmax)
    end
    return 0
end

function contourflocal(f::Field{ProductSpace{S1, S2}}, levels, pmin, pmax) where {S1, S2} 
    u  = Field(f.space.S1, u->u)
    v  = Field(f.space.S2, v->v)
    contourf(v.value, u.value, f.value, levels, vmin=pmin, vmax=pmax)
    return 0
end

function Base. maximum(grid::Array{Field})
    localmax = []
    for index in CartesianIndices(grid)
        append!(localmax, maximum(grid[index]))
    end
    return maximum(localmax)
end

function Base. minimum(grid::Array{Field})
    localmin = []
    for index in CartesianIndices(grid)
        append!(localmin, minimum(grid[index]))
    end
    return minimum(localmin)
end

function levels(grid::Array{Field}, nlevels)
    return collect(range(minimum(grid), stop=maximum(grid), length=nlevels))
end

function extractfield(grid::Array{NTuple{3, Field}}, nfield::Int)::Array{Field}
    fieldarray = Array{Field}(undef, size(grid))
    for index in CartesianIndices(fieldarray)
        fieldarray[index] = grid[index][nfield] 
    end
    return fieldarray
end

function extractoutgoingVboundary(A::NTuple{3, Field})::NTuple{3, Field}
    return (Field(A[1].space.S1, A[1].value[:, 1]), 
            Field(A[2].space.S1, A[2].value[:, 1]),
            Field(A[3].space.S1, A[3].value[:, 1]))
end

function extractoutgoingUboundary(A::NTuple{3, Field})::NTuple{3, Field}
    return (Field(A[1].space.S2, A[1].value[1, :]), 
            Field(A[2].space.S2, A[2].value[1, :]),
            Field(A[3].space.S2, A[3].value[1, :]))
end

function computeUinitialdata(A::NTuple{3, Field})::NTuple{3, Field} 
    (a, r, ϕ) = A
    rU = solveR(extractUboundary(a), extractUboundary(r), extractUboundary(ϕ))
    return (extractUboundary(a), rU, extractUboundary(ϕ))
end

function computeVinitialdata(A::NTuple{3, Field})::NTuple{3, Field}
    (a, r, ϕ) = A
    rV = solveR(extractVboundary(a), extractVboundary(r), extractVboundary(ϕ))
    return (extractVboundary(a), rV, extractVboundary(ϕ))
end

function combineUVboundarylocal(uboundary::NTuple{3, Field{S2}},
                                vboundary::NTuple{3, Field{S1}})::NTuple{3, Field{ProductSpace{S1, S2}}} where {S1, S2}
    aU, rU, ϕU = uboundary
    aV, rV, ϕV = vboundary
    return (combineUVboundary(aU, aV), 
            combineUVboundary(rU, rV), 
            combineUVboundary(ϕU, ϕV))
end

function distribute(ubounds::NTuple{2, Number}, vbounds::NTuple{2, Number}, npatches::Int, npoints::Int)
    ustops = range(ubounds[1], stop=ubounds[2], length=npatches+1) 
    vstops = range(vbounds[1], stop=vbounds[2], length=npatches+1) 
    grid   = Array{NTuple{3, Field}}(undef, npatches, npatches)

    for index in CartesianIndices(grid)
        PS = ProductSpace(ChebyshevGL{U, npoints, Float64}(ustops[index.I[1]], ustops[index.I[1]+1]),
                          ChebyshevGL{V, npoints, Float64}(vstops[index.I[2]], vstops[index.I[2]+1]))
        grid[index] = initialguess(PS, (u,v)->0.6*(exp(-(v-1)^2)))
    end

    for index in CartesianIndices(grid)
        @show index
        uboundary = index.I[1] == 1 ?  computeUinitialdata(grid[index]) : extractoutgoingUboundary(grid[index])
        vboundary = index.I[2] == 1 ?  computeVinitialdata(grid[index]) : extractoutgoingVboundary(grid[index])
        grid[index]  = nonlinearsolver(grid[index][1].space, combineUVboundarylocal(uboundary, vboundary), grid[index])
    end
        
    return grid
end

grid = distribute((0,1), (2,3), 2, 20)
contourfgrid(extractfield(grid, 3), 10)
show()

