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
        return (mix!(mix!(F1, A, FonAxis1), B, a-bnda), 
                mix!(mix!(F2, A, FonAxis2), B, η-bndη))
    end

    function f!(fvec::Array{T,1}, x::Array{T,1}) where {T}
        fvec[:] = reshapeFromTuple(F(reshapeToTuple(PS, x)...))
    end

    
    function j!(jvec::Array{T,2}, x::Array{T,1}) where {T}
        jvec[:, :] = ForwardDiff.jacobian(f!, similar(x), x)
    end

    (bnda, bndη) = boundarydata
    (a0, η0) = initialguess

    solved = reshapeToTuple(PS, nlsolve(f!, reshapeFromTuple((a0, η0)); 
                                        method=:trust_region, factor=0.001, autodiff=:forward,
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
@show L2(asol - a0)
@show L2(ηsol - η0)

