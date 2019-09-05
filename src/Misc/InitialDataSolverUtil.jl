#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 09-2019
# Utilities for initial data solver
#--------------------------------------------------------------------

export initialdatasolver

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

function solveR(f::Field{S}, r::Field{S}, ϕ::Field{S})::Field{S} where {S}
    D = derivative(ϕ.space)
    I = identity(ϕ.space)
    B = incomingboundary(ϕ.space) + outgoingboundary(ϕ.space)
    A = 2*(D*D - (1/f)*(D*f)*(D)) + ((D*ϕ)^2)*I
    return solve(A ⊕ B, B*r)
end

function initialdatasolver(f::Field{S}, r::Field{S}, ϕ::Field{S})::NTuple{3, Field{S}} where {S}
    rU = solveR(extractUboundary(f), extractUboundary(r), extractUboundary(ϕ))
    rV = solveR(extractVboundary(f), extractVboundary(r), extractVboundary(ϕ))
    r  = combineUVboundary(rU, rV)
    return (f, r, ϕ)
end

