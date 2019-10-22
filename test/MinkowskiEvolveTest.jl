#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2019
# Simulate Minkowski spacetime on axis
#--------------------------------------------------------------------

using NLsolve, ForwardDiff
import LinearAlgebra.eigvals, LinearAlgebra.display, LinearAlgebra.cond

function reshapeFromTuple(U::NTuple{2, Field})
    return vcat(reshape(U[1]), reshape(U[2]))
end

function reshapeToTuple(space::S, x::Array{T,1})::NTuple{2, Field}  where {S, T}
    U = reshape(x, :, 2)
    return (reshape(space, U[:, 1]), reshape(space, U[:, 2]))
end

function minkowskisolver(boundarydata::NTuple{2, Field{ProductSpace{S1, S2}}}, 
                         initialguess::NTuple{2, Field{ProductSpace{S1, S2}}})::NTuple{2, Field{ProductSpace{S1, S2}}} where {S1, S2}

    function F(a::Field{S}, η::Field{S})::NTuple{2, Field{S}} where {S}
        F1 = DU*(DV*a) - (1/a)*(DU*a)*(DV*a) + (a/η)*(DU*(DV*η))
        F2 = DU*(DV*η) + (1/η)*(DU*η)*(DV*η) + (1/4)*(1/η)*(a^2)
        FonAxis1 = DU*DV*a - (1/a)*(DU*a)*(DV*a)
        FonAxis2 = DU*DV*η 
        # FonAxis1 = DU*a + DV*a
        # FonAxis2 = DU*η + DV*η 
        # return ((I-B)*(replaceNaNs(F1) + A*FonAxis1) + B*(a-bnda),
                # (I-B)*(replaceNaNs(F2) + A*FonAxis2) + B*(η-bndη))
        return (mix!(mix!(F1, A, FonAxis1), B, a-bnda), 
                mix!(mix!(F2, A, FonAxis2), B, η-bndη))
    end

    # function J(a::Field{S}, η::Field{S})::NTuple{4, Operator{S}} where {S}
        # L1Δa = DU*DV - (1/a)*(DV*a)*DU - (1/a)*(DU*a)*DV + (1/a^2)*(DU*a)*(DV*a)*I + (1/η)*(DU*(DV*η))*I
        # L2Δa = (a/2)*(1/η)*I
        # L1Δη = (a/η)*DU*DV - (a/η^2)*(DU*(DV*η))*I
        # L2Δη = DU*DV + (1/η)*(DV*η)*DU + (1/η)*(DU*η)*DV - (1/η^2)*(DU*η)*(DV*η)*I - (1/4)*(a^2)*(1/η^2)*I
        # D1   = DU*DV - (1/a)*(DV*a)*DU - (1/a)*(DU*a)*DV + (1/a^2)*(DU*a)*(DV*a)*I
        # D2   = DV * DU
        # return ((I-B)*(replaceNaNs(L1Δa) +   A*D2) + B, (I-B)*(replaceNaNs(L1Δη) + 0*A*D1),
                # (I-B)*(replaceNaNs(L2Δa) + 0*A*D1),     (I-B)*(replaceNaNs(L2Δη) +   A*D2) + B)
    # end

    function f!(fvec::Array{T,1}, x::Array{T,1}) where {T}
        fvec[:] = reshapeFromTuple(F(reshapeToTuple(PS, x)...))
    end

    
    function j!(jvec::Array{T,2}, x::Array{T,1}) where {T}
        # jvec[:, :] = reshapeFromTuple2E(J(reshapeToTuple2E(PS, x)...))
        jvec[:, :] = ForwardDiff.jacobian(f!, similar(x), x)
        @show cond(jvec)
        # @show eigvals(jvec)
        # display(jvec)
        # println()
        # println()
    end

    (bnda, bndη) = boundarydata
    (a0, η0) = initialguess

    solved = reshapeToTuple(PS, nlsolve(f!, j!, reshapeFromTuple((a0, η0)); method=:trust_region, factor=0.001,
                                            show_trace=true, ftol=1e-9, iterations=100).zero)

    return solved

end

PS = ProductSpace(ChebyshevGL{U, 8, Float64}(0, 1), 
                  ChebyshevGL{V, 8, Float64}(0, 1))

DU, DV = derivative(PS)
B = incomingboundary(PS)
A = axisboundary(PS)
I = identity(PS)

a0 = Field(PS, (u,v)->1)
η0 = Field(PS, (u,v)->(v-u)/2)
(asol, ηsol) = minkowskisolver((B*a0, B*η0), (0.1 + a0, η0))
display(asol)
display(asol - a0)
@show L2(asol - a0)
@show L2(ηsol - η0)
pcolormesh(asol)
show()
pcolormesh(ηsol)
show()

