#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2019
# Simulate Minkowski spacetime on axis
#--------------------------------------------------------------------

function axisboundary(PS::ProductSpace{S1, S2})::Operator{ProductSpace{S1, S2}} where {S1, S2}
    @assert size(PS.S1) == size(PS.S2) 
    @assert range(PS.S1) == range(PS.S2) 
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

PS = ProductSpace(ChebyshevGL{U, 2, Float64}(0, 1), 
                  ChebyshevGL{V, 2, Float64}(0, 1))

A = axisboundary(PS)
I = identity(PS)

r = Field(PS, (u,v)->2/(v-u))
w = Field(PS, (u,v)->4)

display(r)
display(replaceNaNs(r) + A*w)

