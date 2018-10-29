#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
# TODO: Define the space over which the fields are computed
# TODO: Add assert statements to check if the coordinate boundaries 
#       are satisfied
#--------------------------------------------------------------------

struct Patch
    space::Type{S}
    value::Field
end

struct Database{D1, D2}
    patch::Array{Union{Nothing, Patch}}(nothing, D1, D2)
end

function Database(D::Int) = Database(D, D) 

function neighbour(loc, LR::Symbol)::loc
    if LR==:L
        return loc - [0, 1]
    elseif LR == :R
        return loc - [1, 0]
    else
        @error "Invalid neighbour specified"
    end
end

function getPatchIC(map::Function, u::Field{S}, v::Field{S}, LR::Symbol)::Boundary where {S}
    return Boundary(boundary(Null, SUV, LR)*Field(u.space, map, u, v))
end

function getPatchBnd(database::Database, loc::Array{Int,1}, LR::Symbol)::Boundary
    locboundary = neighbour(loc, LR)
    return Boundary(boundary(Null, :outgoing, LR)*database.patch[locboundary])
end

function distribute(LBoundary::Function, RBoubdary::Function, P1::T, P2::T, M::T)::Database where {T<:Int}
    database = Database(M) 
    guu = Field(SUV, (u,v)-> 0)
    guv = Field(SUV, (u,v)->-2)
    gvv = Field(SUV, (u,v)-> 0)
        
    for i in 2:2M, k in i-min(i-1,M):min(i-1,M) 
        loc = [k, i-k]
        SUV = ProductSpace{GaussLobatto{U,P1}, GaussLobatto{V,P2}}
        ul  = Field(SUV, (u,v)->u)
        vl  = Field(SUV, (u,v)->v)
        Ω   = Field(SUV, (ul,vl)->(pi/8)*cospi(ul/2)*cospi(vl/2))
        
        𝒖 =  ul*cos(Ω) + vl*sin(Ω)
        𝒗 = -ul*sin(Ω) + vl*cos(Ω)
        𝔻𝒗, 𝔻𝒖 = derivativetransform(SUV, 𝒖, 𝒗)
        𝔹      = boundary(Null, SUV)
        
        (𝕘𝒖𝒖, 𝕘𝒖𝒗, 𝕘𝒗𝒗) = inversemetrictransform(guu, guv, gvv, 𝒖, 𝒗) 
        𝕃   = 𝕘𝒖𝒗*𝔻𝒖*𝔻𝒗 + 𝕘𝒖𝒗*𝔻𝒗*𝔻𝒖

        # find the boundaries from neighbouring patches 
        boundaryL = loc[2]==1 ? getPatchIC(LB, 𝒖, 𝒗, :L) : getPatchBnd(database, loc, :L)
        boundaryR = loc[1]==1 ? getPatchIC(RB, 𝒖, 𝒗, :R) : getPatchBnd(database, loc, :R)
        𝕓 = boundaryL + boundaryR

        𝕨 = solve(𝕃 + 𝔹, ρ + 𝕓) 
        database[loc] = Patch(SUV, 𝕨)
    end

    return database
end

