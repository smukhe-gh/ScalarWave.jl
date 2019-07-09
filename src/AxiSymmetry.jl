#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 04-2019
# Define special operators for axis-symmetry
#--------------------------------------------------------------------

export axisboundary, enforceregularityonaxis

function axisboundary(PS::ProductSpace{S1, S2})::Operator{ProductSpace{S1,  S2}} where {S1, S2}
    @assert size(PS.S1) == size(PS.S2)
    return reshape(PS, Diagonal(vec(reshape(identity(PS.S1))))) 
end

function enforceregularityonaxis(A::Operator{ProductSpace{S1, S2}}, D::Operator{ProductSpace{S1, S2}})::Operator{ProductSpace{S1, S2}} where {S1, S2}
    @assert range(A.space) == range(D.space)
    A = axisboundary(A.space)
    return ⊕(A,B,D)
end

