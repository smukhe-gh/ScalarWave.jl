#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2019
# Test analytic and AD jacobian
#--------------------------------------------------------------------

function J(a::Field{S}, η::Field{S})::NTuple{4, Operator{S}} where {S}
    L1Δη = DU*DV + (1/η)*(DV*η)*DU + (1/η)*(DU*η)*DV - (1/η^2)*(DU*η)*(DV*η)*I - (1/4)*(a^2)*(1/η^2)*I
    L1Δa = (a/2)*(1/η)*I
    L2Δη = (a/η)*DU*DV - (a/η^2)*(DU*(DV*η))*I
    L2Δa = DU*DV - (1/a)*(DV*a)*DU - (1/a)*(DU*a)*DV + (1/a^2)*(DU*a)*(DV*a)*I + (1/η)*(DU*(DV*η))*I
    D    = DV + DU
    return ((I-B)*replaceNaNs(L1Δa) + B, (I-B)*(replaceNaNs(L1Δη) + A*D),
            (I-B)*replaceNaNs(L2Δa),     (I-B)*(replaceNaNs(L2Δη) + A*D) + B)
end

function F(a::Field{S}, η::Field{S})::NTuple{2, Field{S}} where {S}
    F1 = DU*(DV*η) + (1/η)*(DU*η)*(DV*η) + (1/4)*(1/η)*(a^2)
    F2 = DU*(DV*a) - (1/a)*(DU*a)*(DV*a) + (a/η)*(DU*(DV*η))
    FonAxis = DV*η + DU*η 
    return ((I-B)*(replaceNaNs(F1) + A*FonAxis) + B*(a-bnda),
            (I-B)*(replaceNaNs(F2) + A*FonAxis) + B*(η-bndη))
end

function f!(fvec::Array{T,1}, x::Array{T,1}) where {T}
    fvec[:] = reshapeFromTuple2E(F(reshapeToTuple2E(PS, x)...))
end

function j!(jvec::Array{T,2}, x::Array{T,1}) where {T}
    jvec[:, :] = reshapeFromTuple2E(J(reshapeToTuple2E(PS, x)...))
end

#--------------------------------------------------------------------
# Test Jacobians
#--------------------------------------------------------------------
PS = ProductSpace(ChebyshevGL{U, 13, Float64}(0, 1), 
                  ChebyshevGL{V, 13, Float64}(0, 1))

DU, DV = derivative(PS)
B = incomingboundary(PS)
I = identity(PS)
A = axisboundary(PS)

a = Field(PS, (u,v)->1)
η = Field(PS, (u,v)->(v-u)/2)
bnda = B*a
bndη = B*η

using ForwardDiff
JfromAD = ForwardDiff.jacobian(f!, reshapeFromTuple2E((a, η)), reshapeFromTuple2E((a, η)))
Janalytic = zeros(2 .*size(PS) .*size(PS))
j!(Janalytic, reshapeFromTuple2E((a, η)))

# println("Automatic differentiation")
# display(JfromAD)
# println()
# println("Analytic")
# display(Janalytic)
# display(JfromAD .== Janalytic)

@test JfromAD ≈ Janalytic
exit()
