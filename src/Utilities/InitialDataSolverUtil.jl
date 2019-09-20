#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 09-2019
# Utilities for initial data solver
#--------------------------------------------------------------------

export initialdatasolver
export extractUboundary, extractVboundary, combineUVboundary, solveR, initialdatasolver, initialguess

function extractUboundary(u::Field{ProductSpace{S1,S2}})::Field{S2} where {S1, S2}
    @assert ndims(u.value) == 2
    return Field(u.space.S2, u.value[end, :])
end

function extractVboundary(u::Field{ProductSpace{S1,S2}})::Field{S1} where {S1, S2}
    @assert ndims(u.value) == 2
    return Field(u.space.S1, u.value[:, end])
end

function combineUVboundary(ubnd::Field{SV}, vbnd::Field{SU})::Field{ProductSpace{SU, SV}} where {SU, SV}
    PS = ProductSpace(vbnd.space, ubnd.space)
    w  = Field(PS, (u,v)->rand())    
    w.value[end, :] = ubnd.value
    w.value[:, end] = vbnd.value
    return w
end

function solveR(a::Field{S}, r::Field{S}, ϕ::Field{S})::Field{S} where {S<:Space{Tag}} where {Tag}
    D = derivative(ϕ.space)
    I = identity(ϕ.space)
    B = incomingboundary(ϕ.space) + outgoingboundary(ϕ.space)
    A = D*D - (2/a)*(D*a)*D + 4pi*(D*ϕ)^2*I
    return solve(A ⊕ B, B*r)
end

function initialdatasolver(a::Field{S}, r::Field{S}, ϕ::Field{S})::NTuple{3, Field{S}} where {S}
    rU = solveR(extractUboundary(a), extractUboundary(r), extractUboundary(ϕ))
    rV = solveR(extractVboundary(a), extractVboundary(r), extractVboundary(ϕ))
    rs  = combineUVboundary(rU, rV)
    B  = incomingboundary(a.space)
    return (B*a, B*rs, B*ϕ)
end

function initialguess(PS::ProductSpace{S1, S2}, map::Function)::NTuple{3, Field{ProductSpace{S1, S2}}} where {S1, S2}
    ϕ0 = Field(PS, map)
    r0 = Field(PS, (u,v)->v-u)
    f0 = Field(PS, (u,v)->1)
    a0 = -sqrt(2*f0) 
    return (a0, r0, ϕ0)
end

