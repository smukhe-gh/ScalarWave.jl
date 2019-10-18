#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2019
# Test analytic and AD jacobian
#--------------------------------------------------------------------

function axisboundary(PS::ProductSpace{S1, S2})::Operator{ProductSpace{S1, S2}} where {S1, S2}
    C = 0*identity(PS)
    for index in CartesianIndices(C.value)
        if index.I[1] == index.I[2] == index.I[3] == index.I[4] 
            C.value[index] = 1 
        end
    end
    return C
end

function replaceNaNs(u::Field{S})::Field{S} where {S}
    return Field(u.space, isfinite.(u.value).*u.value)
end

function replaceNaNs(A::Operator{S})::Operator{S} where {S}
    return Operator(A.space, isfinite.(A.value).*A.value)
end

function J(a::Field{S}, η::Field{S})::NTuple{4, Operator{S}} where {S}
    L1Δη = η*DU*DV + (DV*η)*DU + (DU*η)*DV + (DU*(DV*η))*I
    L1Δa = (a/2)*I
    L2Δη = (1/η)*DU*DV - (1/η^2)*(DU*(DV*η))*I
    L2Δa = (1/a)*DU*DV - (1/a^2)*(DV*a)*DU - (1/a^2)*(DU*a)*DV - (1/a^2)*(DU*(DV*a))*I + (2/a^3)*(DU*a)*(DV*a)*I
    D    = DV + DU
    return ((I-B)*L1Δa + B,    (I-B)*(replaceNaNs(L1Δη) + A*D),
            (I-B)*(I-A)*L2Δa,  (I-B)*(replaceNaNs(L2Δη) + A*D) + B)
end

function F(a::Field{S}, η::Field{S})::NTuple{2, Field{S}} where {S}
    F1 = η*(DU*(DV*η)) + (DU*η)*(DV*η) + (1/4)*(a^2)
    F2 = (1/a)*(DU*(DV*a)) - (1/a^2)*(DU*a)*(DV*a) + (1/η)*(DU*(DV*η))
    FonAxis = DV*η + DU*η 
    return ((I-B)*(replaceNaNs(F1) + A*FonAxis) + B*(a-bnda),
            (I-B)*(replaceNaNs(F2) + A*FonAxis) + B*(η-bndη))
end

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
    A = [reshape(U[1]) reshape(U[2]); 
         reshape(U[3]) reshape(U[4])] 
    return A
end

#--------------------------------------------------------------------
# Test Jacobians
#--------------------------------------------------------------------
PS = ProductSpace(ChebyshevGL{U, 17, Float64}(0, 1), 
                  ChebyshevGL{V, 17, Float64}(0, 1))

DU, DV = derivative(PS)
B = incomingboundary(PS)
I = identity(PS)
A = axisboundary(PS)

a = Field(PS, (u,v)->1)
η = Field(PS, (u,v)->(v-u)/2)
bnda = B*a
bndη = B*η

using ForwardDiff
JfromAD = ForwardDiff.jacobian(f!, reshapeFromTuple((a, η)), reshapeFromTuple((a, η)))
Janalytic = zeros(2 .*size(PS) .*size(PS))
j!(Janalytic, reshapeFromTuple((a, η)))

# println("Automatic differentiation")
# display(JfromAD)
# println()
# println("Analytic")
# display(Janalytic)
# display(JfromAD .== Janalytic)
@test JfromAD ≈ Janalytic
@show cond(JfromAD)
@show cond(Janalytic)
exit()
