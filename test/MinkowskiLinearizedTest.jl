#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2019
# Linearize Einstien's equations on a Minkowski background
#--------------------------------------------------------------------

function regularize(A::Operator{ProductSpace{S1, S2}}, B::Operator{ProductSpace{S1, S2}})::Operator{ProductSpace{S1, S2}} where {S1, S2}
    @assert size(A.space.S1) == size(A.space.S2) 
    @assert range(A.space.S1) == range(A.space.S2) 
    C = Operator(A.space, deepcopy(A.value))
    for index in CartesianIndices(C.value)
        if index.I[1] == index.I[2] 
            C.value[index] = B.value[index]
        end
    end
    return C
end

function L1Δη(PS::ProductSpace{S1, S2})::Operator{ProductSpace{S1, S2}} where {S1, S2}
    return η*DU*DV + (DV*η)*DU + (DU*η)*DV + (DU*(DV*η))*I
end

function L1Δa(PS::ProductSpace{S1, S2})::Operator{ProductSpace{S1, S2}} where {S1, S2}
    return (a/2)*I
end

function L2Δη(PS::ProductSpace{S1, S2})::Operator{ProductSpace{S1, S2}} where {S1, S2}
    return (1/η)*DU*DV - (1/η^2)*(DU*(DV*η))*I
end

function L2Δa(PS::ProductSpace{S1, S2})::Operator{ProductSpace{S1, S2}} where {S1, S2}
    return (1/a)*DU*DV - (1/a^2)*(DV*a)*DU - (1/a^2)*(DU*a)*DV - (1/a^2)*(DU*(DV*a))*I + (2/a^3)*(DU*a)*(DV*a)*I
end

function linearizedEinstein(PS::ProductSpace{S1, S2})::Array{Number, 2} where {S1, S2}
    return [reshape(regularize(L1Δη(PS), D)) reshape(L1Δa(PS)); 
            reshape(regularize(L2Δη(PS), D)) reshape(L2Δa(PS))] 
end

function boundaryconditions(PS::ProductSpace{S1, S2})::Array{Number, 2} where {S1, S2}
    return [reshape(B)   reshape(0*I); 
            reshape(0*I)   reshape(B)] 
end

PS = ProductSpace(ChebyshevGL{U, 10, Float64}(0, 1), 
                  ChebyshevGL{V, 10, Float64}(0, 1))

DU, DV = derivative(PS)
I = identity(PS)
B = incomingboundary(PS)
D = DV + DU
Q = DV - DU

ω = Field(PS, (u,v)->1)
r = Field(PS, (u,v)->(v-u)/2)
a = Field(PS, (u,v)->1)
η = r*ω


A = regularize(L1Δη(PS), D) ⊕ B
@show cond(A)

C = regularize(L2Δη(PS), D) ⊕ B
@show cond(C)

E = L1Δa(PS) ⊕ B
@show cond(E)

F = L2Δa(PS) ⊕ B
@show cond(F)

G = linearizedEinstein(PS)
H = boundaryconditions(PS)
J = G + H

import LinearAlgebra.display, LinearAlgebra.eigvals, LinearAlgebra.rank

@show eigvals(J)
@show cond(J)
@show size(J)
@show rank(J)
