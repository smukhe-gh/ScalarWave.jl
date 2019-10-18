#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 04-2019
# Define operators for axis-symmetry
#--------------------------------------------------------------------

export axisboundary, replaceNaNs

function axisboundary(PS::ProductSpace{S1, S2})::Operator{ProductSpace{S1, S2}} where {S1, S2}
    C = 0*identity(PS)
    for index in CartesianIndices(C.value)
        if index.I[1] == index.I[2] == index.I[3] == index.I[4] 
            C.value[index] = 1 
        end
    end
    return C
end

function replaceNaNs(u::Field{S})::Field{S} where {S}
    return Field(u.space, isfinite.(u.value).*u.value)
end

function replaceNaNs(A::Operator{S})::Operator{S} where {S}
    return Operator(A.space, isfinite.(A.value).*A.value)
end

