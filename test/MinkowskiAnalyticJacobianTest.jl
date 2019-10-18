#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2019
#--------------------------------------------------------------------

using NLsolve

function regularizeJacobian(A::Operator{ProductSpace{S1, S2}}, B::Operator{ProductSpace{S1, S2}})::Operator{ProductSpace{S1, S2}} where {S1, S2}
    @assert size(A.space.S1) == size(A.space.S2) 
    @assert range(A.space.S1) == range(A.space.S2) 
    C = Operator(A.space, deepcopy(A.value))
    for index in CartesianIndices(C.value)
        if index.I[1] == index.I[2] == index.I[3] == index.I[4]
            C.value[index] = B.value[index]
        end
    end
    return C
end

function minkowskisolver(boundarydata::NTuple{2, Field{ProductSpace{S1, S2}}}, 
                         initialguess::NTuple{2, Field{ProductSpace{S1, S2}}})::NTuple{2, Field{ProductSpace{S1, S2}}} where {S1, S2}

    function J(a::Field{S}, η::Field{S})::NTuple{4, Operator{S}} where {S}
        L1Δη = η*DU*DV + (DV*η)*DU + (DU*η)*DV + (DU*(DV*η))*I
        L1Δa = (a/2)*I
        L2Δη = (1/η)*DU*DV - (1/η^2)*(DU*(DV*η))*I
        L2Δa = (1/a)*DU*DV - (1/a^2)*(DV*a)*DU - (1/a^2)*(DU*a)*DV - (1/a^2)*(DU*(DV*a))*I + (2/a^3)*(DU*a)*(DV*a)*I
        D    = DV + DU
        return ((I-B)*regularizeJacobian(L1Δη, D) + B, (I-B)*L1Δa, 
                (I-B)*regularizeJacobian(L2Δη, D)    , (I-B)*L2Δa + B)
    end

    function F(a::Field{S}, η::Field{S})::NTuple{2, Field{S}} where {S}

        function regularizeResidual(F::Field{S})::Field{S} where {S}
            FonAxis = DV*η + DU*η 
            Fcopy = Field(F.space, deepcopy(F.value))
            for index in CartesianIndices(Fcopy.value)
                if index.I[1] == index.I[2] && η.value[index] <= eps(Float64)
                    Fcopy.value[index] = FonAxis.value[index]
                end
            end
            return Fcopy
        end

        F1 = η*(DU*(DV*η)) + (DU*η)*(DV*η) + (1/4)*(a^2)
        F2 = (1/a)*(DU*(DV*a)) - (1/a^2)*(DU*a)*(DV*a) + (1/η)*(DU*(DV*η))

        return (regularizeResidual(F1), regularizeResidual(F2))
    end

    function C(a::Field{S}, η::Field{S})::NTuple{2, Field{S}} where {S}
        DU, DV = derivative(a.space)
        C1 = DU*(DU*η) - (2/a)*(DU*a)*(DU*η) 
        C2 = DV*(DV*η) - (2/a)*(DV*a)*(DV*η) 
        return (C1, C2)
    end

    function residual(a::Field{S}, η::Field{S})::NTuple{2, Field{S}} where {S}
        function BCs(F1::Field{S}, F2::Field{S})::NTuple{2, Field{S}} where {S}
            BF1 = (I-B)*F1 + B*(a-bnda) 
            BF2 = (I-B)*F2 + B*(η-bndη) 
            return (BF1, BF2)
        end
        return BCs(F(a, η)...)
    end


    function f!(fvec::Array{T,1}, x::Array{T,1}) where {T}
        fvec[:] = reshapeFromTuple(residual(reshapeToTuple(PS, x)...))
    end

    function j!(jvec::Array{T,2}, x::Array{T,1}) where {T}
        jvec[:, :] = reshapeFromTuple(J(reshapeToTuple(PS, x)...))
    end

    function reshapeFromTuple(U::NTuple{2, Field})
        return vcat(reshape(U[1]), reshape(U[2]))
    end
    
    function reshapeToTuple(space::S, x::Array{T,1})::NTuple{2, Field}  where {S, T}
        U = reshape(x, :, 2)
        return (reshape(space, U[:, 1]), reshape(space, U[:, 2]))
    end

    function reshapeFromTuple(U::NTuple{4, Operator})
        return [reshape(U[1]) reshape(U[2]); 
                reshape(U[3]) reshape(U[4])] 
    end

    (bnda, bndη) = boundarydata
    (a0, η0) = initialguess

    @show L2.(C(initialguess...))
    @show L2.(F(initialguess...))

    solved = reshapeToTuple(PS, nlsolve(f!, j!, reshapeFromTuple((a0, η0)); 
                                            show_trace=true, ftol=1e-9, iterations=50).zero)
    @show L2.(C(solved...))
    @show L2.(F(solved...))
    return solved
end

PS = ProductSpace(ChebyshevGL{U, 16, Float64}(0, 1), 
                  ChebyshevGL{V, 16, Float64}(0, 1))

DU, DV = derivative(PS)
B = incomingboundary(PS)
I = identity(PS)

a0 = Field(PS, (u,v)->1)
η0 = Field(PS, (u,v)->(v-u)/2)

(asol, rsol) = minkowskisolver((B*a0, B*η0), (1 + a0, η0))
pcolormesh(asol)
show()
pcolormesh(rsol)
show()

