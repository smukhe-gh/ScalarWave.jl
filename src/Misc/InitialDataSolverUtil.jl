#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 09-2019
# Utilities for initial data solver
#--------------------------------------------------------------------

export initialdatasolver
export extractUboundary, extractVboundary, combineUVboundary

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
