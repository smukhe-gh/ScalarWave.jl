#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 09-2019
# Utilities for initial data solver
#--------------------------------------------------------------------

export Grid, setup, extractfield
export computeUboundary, computeVboundary, compute, distribute

struct Grid
    eltype::Type
    npoints::NTuple{2, Int}
    npatches::NTuple{2, Int}
    ubounds::NTuple{2, Number}
    vbounds::NTuple{2, Number}
end

function extractfield(tree::Array{NTuple{3, Field}, 2}, ID::Symbol)::Array{Field, 2}
    ID = (ID == :a ? 1 : ( ID == :r ? 2 : 3))
    branch = Array{Field}(undef, size(tree))
    for index in CartesianIndices(branch)
        branch[index] = tree[index][ID] 
    end
    return branch
end

function compute(nonlinearsolver::Function,
                 ubnd::NTuple{3, Field{S2}}, 
                 vbnd::NTuple{3, Field{S1}}, 
                 background::NTuple{3, Field{ProductSpace{S1, S2}}})::NTuple{3, Field{ProductSpace{S1, S2}}} where {S1, S2}
    return nonlinearsolver(combineUVboundary(ubnd, vbnd, :incoming), background)
end
    

function setup(grid::Grid, background::Function)::Array{NTuple{3, Field}}
    ustops = range(grid.ubounds[1], stop=grid.ubounds[2], length=grid.npatches[1]+1) 
    vstops = range(grid.vbounds[1], stop=grid.vbounds[2], length=grid.npatches[2]+1) 
    tree   = Array{NTuple{3, Field}}(undef, grid.npatches[1], grid.npatches[2])
    for index in CartesianIndices(tree)
        PS = ProductSpace(ChebyshevGL{U, grid.npoints[1], grid.eltype}(ustops[index.I[1]], ustops[index.I[1]+1]),
                          ChebyshevGL{V, grid.npoints[2], grid.eltype}(vstops[index.I[2]], vstops[index.I[2]+1]))
        tree[index] = background(PS)
    end
    return tree
end

function distribute(grid::Grid, background::Function, nonlinearsolver::Function)::Array{NTuple{3, Field}}
    tree = setup(grid, background)
    for index in CartesianIndices(tree)
        uboundary = index.I[1] == 1 ?  computeUboundary(tree[index]) : extractUboundary(tree[index - CartesianIndex((1,0))], :outgoing)
        vboundary = index.I[2] == 1 ?  computeVboundary(tree[index]) : extractVboundary(tree[index - CartesianIndex((0,1))], :outgoing)
        tree[index] = compute(nonlinearsolver, uboundary, vboundary, tree[index])
    end
    return tree
end
