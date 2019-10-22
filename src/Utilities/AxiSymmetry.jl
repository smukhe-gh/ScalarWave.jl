#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 04-2019
# Define operators for axis-symmetry
#--------------------------------------------------------------------

export axisboundary, replaceNaNs
export mix!, mixonaxis, mixonincomingboundary

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

function mix!(u::Field{S}, A::Operator{S}, v::Field{S})::Field{S} where {S}
    for index in CartesianIndices(u.value)
        if A.value[index.I..., index.I...] == 1
            u.value[index] = v.value[index]
        end
    end
    return u
end

function mixonaxis(u::Field{S}, v::Field{S})::Field{S} where {S}
    for index in CartesianIndices(u.value)
        if index.I[1] == index.I[2]
            u.value[index] == v.value[index]
        end
    end
    return u
end

function mixonincomingboundary(u::Field{S}, v::Field{S})::Field{S} where {S}
    indU, indV = size(u.space)
    for index in CartesianIndices(u.value)
        if index.I[1] == indU || index.I[2] == indV
            u.value[index] == v.value[index]
        end
    end
    return u
end


# function mix(B::Operator, int::Field, bnd::Feld)
    # r = similar(int)
    # for index in CartesianIndide(r)
        # if B[index, index] == 0
            # r[idnex] = ind[index]
        # else
            # r[index] = bnd[index['
                              # ]
                           # ;
                           # '
                           # [
                           # ]
                          # ]
        # end
        # r
    # end


# function mix(B::Operator, int::Operator, bnd::Operator)
    # r = similar(int)
    # for indiex in CartesianIndices(r)
        # if B[index] == 0
            # r[index] =  int[index]
        # else
            # r[index]  = bnd[index]
        # end
    # end
    # r
# end



