#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 09-2019
# Utilities for nonlinear solver
#--------------------------------------------------------------------

using NLsolve
export reshapeFromTuple, reshapeToTuple, nonlinearsolver

function reshapeFromTuple(U::NTuple{3, Field})
    return vcat(reshape(U[1]), reshape(U[2]), reshape(U[3]))
end

function reshapeToTuple(space::S, x::Array{T,1})::NTuple{3, Field}  where {S, T}
    U = reshape(x, :, 3)
    return (reshape(space, U[:, 1]), reshape(space, U[:, 2]), reshape(space, U[:, 3]))
end

function nonlinearsolver(PS::ProductSpace{S1, S2}, 
                         boundarydata::NTuple{3, Field{ProductSpace{S1, S2}}}, 
                         initialguess::NTuple{3, Field{ProductSpace{S1, S2}}})::NTuple{3, Field{ProductSpace{S1, S2}}} where {S1, S2}

    function C(a::Field{S}, r::Field{S}, ϕ::Field{S})::NTuple{2, Field{S}} where {S}
        C1 = DU*(DU*r) - (2/a)*(DU*a)*(DU*r) + (4pi*r)*(DU*ϕ)^2
        C2 = DV*(DV*r) - (2/a)*(DV*a)*(DV*r) + (4pi*r)*(DV*ϕ)^2
        return (C1, C2)
    end
    
    function F(a::Field{S}, r::Field{S}, ϕ::Field{S})::NTuple{3, Field{S}} where {S}
        F1 = r*(DU*(DV*ϕ)) + (DU*r)*(DV*ϕ) + (DV*r)*(DU*ϕ)
        F2 = r*(DU*(DV*r)) + (DU*r)*(DV*r) + (1/4)*(a^2)
        F3 = (1/a)*(DU*(DV*a)) - (1/a^2)*(DU*a)*(DV*a) + (1/r)*(DU*(DV*r)) + 4pi*(DU*ϕ)*(DV*ϕ)
        return (F1, F2, F3)
    end

    function perfF(a::Field{S}, r::Field{S}, ϕ::Field{S}, DU::Operator{S1}, DV::Operator{S2})::NTuple{3, Field{S}} where {S, S1, S2} # <--optimized with partial sums
        F1 = r.value.*(DU.value*(ϕ.value*transpose(DV.value))) + (DU.value*r.value).*transpose((DV.value*transpose(ϕ.value))) + (DU.value*ϕ.value).*transpose((DV.value*transpose(r.value)))
        F2 = r.value.*(DU.value*(r.value*transpose(DV.value))) + (DU.value*r.value).*transpose((DV.value*transpose(r.value))) + (1/4)*(a.value.*a.value)
        F3 = ((1/a).value.*(DU.value*(a.value*transpose(DV.value))) 
              - (1/a^2).value.*(DU.value*a.value).*transpose((DV.value*transpose(a.value))) 
              + (1/r).value.*(DU.value*(r.value*transpose(DV.value))) 
              + 4pi*(DU.value*ϕ.value).*transpose((DV.value*transpose(ϕ.value))))
        return (Field(a.space, F1), Field(a.space, F2), Field(a.space, F3))
    end
    
    function residual(a::Field{S}, r::Field{S}, ϕ::Field{S})::NTuple{3, Field{S}} where {S}

        function BCs(F1::Field{S}, F2::Field{S}, F3::Field{S})::NTuple{3, Field{S}} where {S}
            BF1 = (I-B)*F1 + B*(a-bnda) 
            BF2 = (I-B)*F2 + B*(r-bndr) 
            BF3 = (I-B)*F3 + B*(ϕ-bndϕ) 
            return (BF1, BF2, BF3)
        end

        return BCs(perfF(a, r, ϕ, dU, dV)...)
    end

    function f!(fvec::Array{T,1}, x::Array{T,1}) where {T}
        fvec[:] = reshapeFromTuple(residual(reshapeToTuple(PS, x)...))
    end

    DU, DV = derivative(PS)
    dU, dV = derivative(PS.S1), derivative(PS.S2)
    B = incomingboundary(PS)
    I = identity(PS)

    (bnda, bndr, bndϕ) = boundarydata
    (a0, r0, ϕ0) = initialguess

    asolved, rsolved, ϕsolved = reshapeToTuple(PS, nlsolve(f!, reshapeFromTuple((a0, r0, ϕ0)); 
                                                           autodiff=:forward, show_trace=true, ftol=1e-9).zero)
    E1, E2 = C(asolved, rsolved, ϕsolved)
    @show L2(E1) + L2(E2)

    return (asolved, rsolved, ϕsolved)
end
