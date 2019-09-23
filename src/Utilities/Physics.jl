#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 09-2019
# Utilities for initial data solver
#--------------------------------------------------------------------

using NLsolve
export background, solveR, lineconstraint, nonlinearsolver, constraints

function background(PS::ProductSpace{S1, S2}, ϕ::Function)::NTuple{3, Field{ProductSpace{S1, S2}}} where {S1, S2}
    ϕ0 = Field(PS, ϕ)
    f0 = Field(PS, (u,v)->1)
    r0 = Field(PS, (u,v)->v-u)
    a0 = -sqrt(2*f0) 
    return (a0, r0, ϕ0)
end

function background(PS::ProductSpace{S1, S2}, f::Function, r::Function, ϕ::Function)::NTuple{3, Field{ProductSpace{S1, S2}}} where {S1, S2}
    ϕ0 = Field(PS, ϕ)
    f0 = Field(PS, f)
    r0 = Field(PS, r)
    a0 = -sqrt(2*f0) 
    return (a0, r0, ϕ0)
end

function solveR(a::Field{S}, r::Field{S}, ϕ::Field{S})::Field{S} where {S<:Space{Tag}} where {Tag}
    D = derivative(ϕ.space)
    I = identity(ϕ.space)
    B = incomingboundary(ϕ.space) + outgoingboundary(ϕ.space)
    A = D*D - (2/a)*(D*a)*D + 4pi*(D*ϕ)^2*I
    return solve(A ⊕ B, B*r)
end

function lineconstraint(a::Field{S}, r::Field{S}, ϕ::Field{S})::Number where {S<:Space{Tag}} where {Tag}
    D = derivative(ϕ.space)
    cres = D*(D*r) - (2/a)*(D*a)*(D*r) + (4pi*r)*(D*ϕ)^2
    return L2(cres)
end

function constraints(a::Field{S}, r::Field{S}, ϕ::Field{S})::NTuple{2, Field{S}} where {S}
    DU, DV = derivative(a.space)
    C1 = DU*(DU*r) - (2/a)*(DU*a)*(DU*r) + (4pi*r)*(DU*ϕ)^2
    C2 = DV*(DV*r) - (2/a)*(DV*a)*(DV*r) + (4pi*r)*(DV*ϕ)^2
    return (C1, C2)
end

function nonlinearsolver(boundarydata::NTuple{3, Field{ProductSpace{S1, S2}}}, 
                         initialguess::NTuple{3, Field{ProductSpace{S1, S2}}})::NTuple{3, Field{ProductSpace{S1, S2}}} where {S1, S2}

    
    function F(a::Field{S}, r::Field{S}, ϕ::Field{S})::NTuple{3, Field{S}} where {S}
        F1 = r*(DU*(DV*ϕ)) + (DU*r)*(DV*ϕ) + (DV*r)*(DU*ϕ)
        F2 = r*(DU*(DV*r)) + (DU*r)*(DV*r) + (1/4)*(a^2)
        F3 = (1/a)*(DU*(DV*a)) - (1/a^2)*(DU*a)*(DV*a) + (1/r)*(DU*(DV*r)) + 4pi*(DU*ϕ)*(DV*ϕ)
        return (F1, F2, F3)
    end

    function residual(a::Field{S}, r::Field{S}, ϕ::Field{S})::NTuple{3, Field{S}} where {S}

        function BCs(F1::Field{S}, F2::Field{S}, F3::Field{S})::NTuple{3, Field{S}} where {S}
            BF1 = (I-B)*F1 + B*(a-bnda) 
            BF2 = (I-B)*F2 + B*(r-bndr) 
            BF3 = (I-B)*F3 + B*(ϕ-bndϕ) 
            return (BF1, BF2, BF3)
        end

        return BCs(F(a, r, ϕ)...)
    end

    function f!(fvec::Array{T,1}, x::Array{T,1}) where {T}
        fvec[:] = reshapeFromTuple(residual(reshapeToTuple(PS, x)...))
    end

    (bnda, bndr, bndϕ) = boundarydata
    (a0, r0, ϕ0) = initialguess

    PS = a0.space
    DU, DV = derivative(PS)
    B = incomingboundary(PS)
    I = identity(PS)

    asolved, rsolved, ϕsolved = reshapeToTuple(PS, nlsolve(f!, reshapeFromTuple((a0, r0, ϕ0)); 
                                                           autodiff=:forward, show_trace=true, ftol=1e-9).zero)
    return (asolved, rsolved, ϕsolved)
end
