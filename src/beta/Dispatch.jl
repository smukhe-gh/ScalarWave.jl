#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
# TODO: Refine Grid to make it more flexible
#--------------------------------------------------------------------

struct Grid
    M::Int
    P::Int
    U::Tuple
    V::Tuple
end

struct Patch{S,D,T}
    space::Type{S}
    value::Field{S,D,T}
end

Database(D::Int) = Array{Union{Nothing, Patch}}(nothing, D, D)
Grid(M, P) = Grid(M, P, (-1, 1), (-1, 1))

function neighbour(loc, LR::Symbol)
    if LR == :L
        return (loc - [0, 1])
    elseif LR == :R
        return (loc - [1, 0])
    else
        @error "Invalid neighbour specified"
    end
end

function getPatchIC(map::Function, u::Field{S}, v::Field{S}, LR::Symbol) where {S}
    return boundary(S, Null, LR)*Field(u.space, map, u, v)
end

function getPatchIC(S::Type{T}, map::Function, LR::Symbol) where {T<:ProductSpace}
    # TODO: Insert assert statements
    return boundary(S, Null, LR)*Field(S, map)
end

function getPatchBnd(S::Type{T}, database::Array, loc::Array{Int,1}, LR::Symbol) where {T<:ProductSpace}
    locboundary = neighbour(loc, LR)
    @show  boundary(S, Null, :outgoing, LR).space
    @show database[locboundary[1], locboundary[2]].value.space
    return boundary(S, Null, :outgoing, LR)*(database[locboundary[1], locboundary[2]].value)
end

function PatchSpace(loc::Array{Int, 1}, grid::Grid)
    umin = range(grid.U[1], stop=grid.U[2], length=grid.M+1)[loc[1]]
    umax = range(grid.U[1], stop=grid.U[2], length=grid.M+1)[loc[1]+1]
    vmin = range(grid.V[1], stop=grid.V[2], length=grid.M+1)[loc[2]]
    vmax = range(grid.V[1], stop=grid.V[2], length=grid.M+1)[loc[2]+1]
    @assert umin === -1.0
    return ProductSpace{GaussLobatto{V, grid.P, vmax, vmin}, 
                        GaussLobatto{U, grid.P, umax, umin}}
end


function distribute(grid::Grid, LB::Function, 
                                RB::Function)
    database = Database(grid.M) 

    for i in 2:2*grid.M, k in i-min(i-1,grid.M):min(i-1,grid.M) 
        loc = [k, i-k]
        @show loc
        SUV = PatchSpace(loc, grid)
        @show SUV
        boundaryL = loc[2]==1 ? getPatchIC(SUV, LB, :L) : getPatchBnd(SUV, database, loc, :L)
        boundaryR = loc[1]==1 ? getPatchIC(SUV, RB, :R) : getPatchBnd(SUV, database, loc, :R)
        𝕓 = boundaryL + boundaryR

        𝕨 = Field(SUV, (u,v)->k)
        database[loc[1], loc[2]] = Patch(SUV, 𝕨)
    end
    return database
end

