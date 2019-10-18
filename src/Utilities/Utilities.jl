#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 09-2019
# Utilities for initial data solver
#--------------------------------------------------------------------

export extractUboundary, extractVboundary, combineUVboundary
export reshapeFromTuple, reshapeToTuple
export reshapeFromTuple2E, reshapeToTuple2E
export Grid, setup, extractfield
export computeUboundary, computeVboundary, compute, distribute

struct Grid
    eltype::Type
    npoints::NTuple{2, Int}
    npatches::NTuple{2, Int}
    ubounds::NTuple{2, Number}
    vbounds::NTuple{2, Number}
end


function reshapeFromTuple(U::NTuple{3, Field})
    return vcat(reshape(U[1]), reshape(U[2]), reshape(U[3]))
end

function reshapeToTuple(space::S, x::Array{T,1})::NTuple{3, Field}  where {S, T}
    U = reshape(x, :, 3)
    return (reshape(space, U[:, 1]), reshape(space, U[:, 2]), reshape(space, U[:, 3]))
end

function extractUboundary(u::Field{ProductSpace{S1, S2}}, boundarytype::Symbol)::Field{S2} where {S1, S2}
    @assert ndims(u.value) == 2
    return (boundarytype == :incoming ? Field(u.space.S2, u.value[end, :]) : Field(u.space.S2, u.value[1, :]))
end

function extractVboundary(u::Field{ProductSpace{S1, S2}}, boundarytype::Symbol)::Field{S1} where {S1, S2}
    @assert ndims(u.value) == 2
    return (boundarytype == :incoming ? Field(u.space.S1, u.value[:, end]) : Field(u.space.S1, (u.value[:, 1])))
end

function combineUVboundary(uboundary::Field{S2}, vboundary::Field{S1}, boundarytype::Symbol)::Field{ProductSpace{S1, S2}} where {S1, S2}
    PS = ProductSpace(vboundary.space, uboundary.space)
    w = Field(PS, (u,v)->0)    
    if boundarytype == :incoming
        w.value[end, :] = uboundary.value
        w.value[:, end] = vboundary.value
    else 
        w.value[1, :]   = uboundary.value
        w.value[:, 1]   = vboundary.value
    end
    return w
end

function extractUboundary(u::NTuple{3, Field{ProductSpace{S1, S2}}}, boundarytype::Symbol)::NTuple{3, Field{S2}} where {S1, S2}
    return extractUboundary.(u, boundarytype) 
end

function extractVboundary(u::NTuple{3, Field{ProductSpace{S1, S2}}}, boundarytype::Symbol)::NTuple{3, Field{S1}} where {S1, S2}
    return extractVboundary.(u, boundarytype) 
end

function combineUVboundary(ubnd::NTuple{3, Field{S2}},
                           vbnd::NTuple{3, Field{S1}}, boundarytype::Symbol)::NTuple{3, Field{ProductSpace{S1, S2}}} where {S1, S2}
    return combineUVboundary.(ubnd, vbnd, boundarytype) 
end

function initialdata(a::Field{S}, r::Field{S}, ϕ::Field{S})::NTuple{3, Field{S}} where {S <: ProductSpace{S1, S2}} where {S1, S2}
    ru = solveR(extractUboundary(a, :incoming), extractUboundary(r, :incoming), extractUboundary(ϕ, :incoming))
    rv = solveR(extractVboundary(a, :incoming), extractVboundary(r, :incoming), extractVboundary(ϕ, :incoming))
    rs = combineUVboundary(ru, rv)
    return (a, rs, ϕ)
end

function computeUboundary(u::NTuple{3, Field{ProductSpace{S1, S2}}})::NTuple{3, Field{S2}} where {S1, S2}
    rs = solveR(extractUboundary.(u, :incoming)...)
    return (extractUboundary(u[1], :incoming), rs, extractUboundary(u[3], :incoming))
end

function computeVboundary(u::NTuple{3, Field{ProductSpace{S1, S2}}})::NTuple{3, Field{S1}} where {S1, S2}
    rs = solveR(extractVboundary.(u, :incoming)...)
    return (extractVboundary(u[1], :incoming), rs, extractVboundary(u[3], :incoming))
end

function extractfield(tree::Array{NTuple{3, Field}, 2}, ID::Symbol)::Array{Field, 2}
    ID = (ID == :a ? 1 : ( ID == :r ? 2 : 3))
    branch = Array{Field}(undef, size(tree))
    for index in CartesianIndices(branch)
        branch[index] = tree[index][ID] 
    end
    return branch
end

function compute(ubnd::NTuple{3, Field{S2}}, 
                 vbnd::NTuple{3, Field{S1}}, background::NTuple{3, Field{ProductSpace{S1, S2}}})::NTuple{3, Field{ProductSpace{S1, S2}}} where {S1, S2}
    return nonlinearsolver(combineUVboundary(ubnd, vbnd, :incoming), background)
end
    
function setup(grid::Grid, ϕ::Function)::Array{NTuple{3, Field}}
    ustops = range(grid.ubounds[1], stop=grid.ubounds[2], length=grid.npatches[1]+1) 
    vstops = range(grid.vbounds[1], stop=grid.vbounds[2], length=grid.npatches[2]+1) 
    tree   = Array{NTuple{3, Field}}(undef, grid.npatches[1], grid.npatches[2])
    for index in CartesianIndices(tree)
        PS = ProductSpace(ChebyshevGL{U, grid.npoints[1], grid.eltype}(ustops[index.I[1]], ustops[index.I[1]+1]),
                          ChebyshevGL{V, grid.npoints[2], grid.eltype}(vstops[index.I[2]], vstops[index.I[2]+1]))
        tree[index] = background(PS, ϕ)
    end
    return tree
end

function distribute(grid::Grid, ϕ::Function)::Array{NTuple{3, Field}}
    tree = setup(grid, ϕ)
    for index in CartesianIndices(tree)
        uboundary = index.I[1] == 1 ?  computeUboundary(tree[index]) : extractUboundary(tree[index - CartesianIndex((1,0))], :outgoing)
        vboundary = index.I[2] == 1 ?  computeVboundary(tree[index]) : extractVboundary(tree[index - CartesianIndex((0,1))], :outgoing)
        tree[index] = compute(uboundary, vboundary, tree[index])
    end
    return tree
end

function reshapeToTuple2E(space::S, x::Array{T,1})::NTuple{2, Field}  where {S, T}
    U = reshape(x, :, 2)
    return (reshape(space, U[:, 1]), reshape(space, U[:, 2]))
end

function reshapeFromTuple2E(U::NTuple{2, Field})
    return vcat(reshape(U[1]), reshape(U[2]))
end

function reshapeFromTuple2E(U::NTuple{4, Operator})
    A = [reshape(U[1]) reshape(U[2]); 
         reshape(U[3]) reshape(U[4])] 
    return A
end
