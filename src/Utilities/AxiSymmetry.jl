#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 04-2019
# Define operators relevant for axis-symmetry
#--------------------------------------------------------------------

export axisboundary, mix! 

function axisboundary(PS::ProductSpace{S1, S2})::Operator{ProductSpace{S1, S2}} where {S1, S2}
    C = 0*identity(PS)
    for index in CartesianIndices(C.value)
        if index.I[1] == index.I[2] == index.I[3] == index.I[4] 
            C.value[index] = 1 
        end
    end
    return C
end

function mix!(u::Field{S}, A::Operator{S}, v::Field{S})::Field{S} where {S}
    for index in CartesianIndices(u.value)
        if A.value[index.I..., index.I...] == 1
            u.value[index] = v.value[index]
        end
    end
    return u
end

