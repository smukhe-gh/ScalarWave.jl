#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2019
# Compute Jacobian 
#--------------------------------------------------------------------

function f!(fvec::Array{T,1}, x::Array{T,1}) where {T}
    fvec[:] = reshapeFromTuple(F(reshapeToTuple(PS, x)...))
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

function J(a::Field{S}, η::Field{S})::NTuple{4, Operator{S}} where {S}
    L1Δη = η*DU*DV + (DV*η)*DU + (DU*η)*DV + (DU*(DV*η))*I
    L1Δa = (a/2)*I
    L2Δη = (1/η)*DU*DV - (1/η^2)*(DU*(DV*η))*I
    L2Δa = (1/a)*DU*DV - (1/a^2)*(DV*a)*DU - (1/a^2)*(DU*a)*DV - (1/a^2)*(DU*(DV*a))*I + (2/a^3)*(DU*a)*(DV*a)*I
    return (L1Δa, L1Δη, 
            L2Δa, L2Δη)
end

function F(a::Field{S}, η::Field{S})::NTuple{2, Field{S}} where {S}
    F1 = η*(DU*(DV*η)) + (DU*η)*(DV*η) + (1/4)*(a^2)
    F2 = (1/a)*(DU*(DV*a)) - (1/a^2)*(DU*a)*(DV*a) + (1/η)*(DU*(DV*η))
    return (F1, F2)
end

#--------------------------------------------------------------------
# Test Jacobian
#--------------------------------------------------------------------

PS = ProductSpace(ChebyshevGL{U, 2, Float64}(2, 3), 
                  ChebyshevGL{V, 2, Float64}(4, 5))

DU, DV = derivative(PS)
B = incomingboundary(PS)
I = identity(PS)

a = Field(PS, (u,v)->4)
η = Field(PS, (u,v)->(v-u)/3)

using ForwardDiff
JfromAD = ForwardDiff.jacobian(f!, reshapeFromTuple((a, η)), reshapeFromTuple((a, η)))

Janalytic = zeros(8,8)
j!(Janalytic, reshapeFromTuple((a, η)))

@test JfromAD ≈ Janalytic
