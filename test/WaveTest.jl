#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2019
# Simulate Minkowski spacetime on axis
#--------------------------------------------------------------------

using NLsolve, ForwardDiff
import LinearAlgebra.eigvals, LinearAlgebra.display, LinearAlgebra.cond

function reshapeFromTuple(U::Field)
    return vcat(reshape(U)) 
end

function reshapeToTuple(space::S, x::Array{T,1})::Field where {S, T}
    U = reshape(x, :, 1)
    return reshape(space, U[:, 1])
end

function einstein(boundarydata::Field{ProductSpace{S1, S2}}, 
                  initialguess::Field{ProductSpace{S1, S2}})::Field{ProductSpace{S1, S2}} where {S1, S2}

    function F(ϕ::Field{S})::Field{S} where {S}
        F3 = DU*(DV*ϕ) + (1/η)*(DU*η)*(DV*ϕ) + (1/η)*(DV*η)*(DU*ϕ) + sin(t) - (1/η)*cos(t)*(DU*η + DV*η) # Add an additional term for zero residual
        F3onAxis = DU*DV*ϕ + sin(t)
        return mix!(mix!(F3, A, F3onAxis), B, ϕ-bndϕ)
    end

    function f!(fvec::Array{T,1}, x::Array{T,1}) where {T}
        fvec[:] = reshapeFromTuple(F(reshapeToTuple(PS, x)))
    end

    
    function j!(jvec::Array{T,2}, x::Array{T,1}) where {T}
        jvec[:, :] = ForwardDiff.jacobian(f!, similar(x), x)
    end

    bndϕ = boundarydata 
    solved = reshapeToTuple(PS, nlsolve(f!, j!, reshapeFromTuple(initialguess); 
                                            show_trace=true, ftol=1e-9, iterations=100).zero)

    return solved
end

PS = ProductSpace(ChebyshevGL{U, 30, Float64}(0, 1), 
                  ChebyshevGL{V, 30, Float64}(0, 1))

DU, DV = derivative(PS)
B = incomingboundary(PS)
A = axisboundary(PS)
I = identity(PS)

t = Field(PS, (u,v)->u+v)
# a = Field(PS, (u,v)->1)
η = Field(PS, (u,v)->(v-u)/2)
ϕ = Field(PS, (u,v)->sin(u+v))
a = 1 + η^2 

ϕsol = einstein(B*ϕ, B*ϕ)
@show L2(ϕsol - ϕ)
