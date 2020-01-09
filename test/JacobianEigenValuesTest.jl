#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2019
# Simulate Minkowski spacetime on axis
#--------------------------------------------------------------------

using NLsolve, ForwardDiff
import LinearAlgebra.eigvals, LinearAlgebra.display, LinearAlgebra.cond

function computeJacobians(PS::ProductSpace{S1, S2}) where {S1, S2}

    function F(a::Field{S}, η::Field{S}, ϕ::Field{S})::NTuple{3, Field{S}} where {S}
        F1 = DU*(DV*a) - (1/a)*(DU*a)*(DV*a) + (a/η)*(DU*(DV*η)) + (4pi*a)*(DU*ϕ)*(DV*ϕ)
        F2 = DU*(DV*η) + (1/η)*(DU*η)*(DV*η) + (1/4)*(1/η)*(a^2)
        F3 = DU*(DV*ϕ) + (1/η)*(DU*η)*(DV*ϕ) + (1/η)*(DV*η)*(DU*ϕ)
        return (mix!(F1, B, a-bnda), 
                mix!(F2, B, η-bndη),
                mix!(F3, B, ϕ-bndϕ))
    end

    function FinRT(a::Field{S}, η::Field{S}, ϕ::Field{S})::NTuple{3, Field{S}} where {S}
        DT =  (1/2)*DU + (1/2)*DV
        DR = -(1/2)*DV + (1/2)*DU

        F1 = (DT*(DT*a) - DR*(DR*a)) - (1/a)*((DT*a)^2 - (DR*a)^2) + (a/η)* (DT*(DT*η) - DR*(DR*η)) + (4pi*a)*(DT*ϕ - DR*ϕ)*(DT*ϕ + DR*ϕ)
        F2 = (DT*(DT*η) - DR*(DR*η)) + (1/η)*((DT*η)^2 - (DR*η)^2) + (1/4)*(1/η)*(a^2)
        F3 = (DT*(DT*ϕ) - DR*(DR*ϕ)) + (1/η)*(DT*η - DR*η)*(DT*ϕ + DR*ϕ) + (1/η)*(DT*η - DR*η)*(DT*ϕ + DR*ϕ)

        return (mix!(F1, B, a-bnda), 
                mix!(F2, B, η-bndη),
                mix!(F3, B, ϕ-bndϕ))
    end

    function f!(fvec::Array{T,1}, x::Array{T,1}) where {T}
        fvec[:] = reshapeFromTuple(FinRT(reshapeToTuple(PS, x)...))
    end

    
    function j!(jvec::Array{T,2}, x::Array{T,1}) where {T}
        jvec[:, :] = ForwardDiff.jacobian(f!, similar(x), x)
    end

    (bnda, bndη, bndϕ) = (B*a, B*η, B*ϕ) 

    J = zeros(3 .* size(a.space).^2)
    x = reshapeFromTuple((a, η, ϕ))
    J = j!(J, x)
    @show sort(abs.(eigvals(J)))
end

PS = ProductSpace(ChebyshevGL{U, 18, Float64}(0, 1), 
                  ChebyshevGL{V, 18, Float64}(0.1, 1.1))

DU, DV = derivative(PS)
B = incomingboundary(PS)
A = axisboundary(PS)
I = identity(PS)

p  = 1e-1
a = Field(PS, (u,v)->1)
η = Field(PS, (u,v)->(v-u)/2)
ϕ = Field(PS, (u,v)->p*exp(-(v-1)^2))

computeJacobians(PS)
