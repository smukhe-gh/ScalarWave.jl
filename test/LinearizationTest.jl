#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2019
# Simulate Minkowski spacetime on axis
#--------------------------------------------------------------------

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
    A = [reshape(U[1]) reshape(U[2]); 
         reshape(U[3]) reshape(U[4])] 
    @show cond(A)
    return A
end

function renormalize(u::Field{S})::Field{S} where {S}
    for index in CartesianIndices(u.value)
        if index.I[1] == index.I[2]
            u.value[index] = 0
        end
    end
    return u
end

function renormalize(A::Operator{S})::Operator{S} where {S}
    for index in CartesianIndices(A.value)
        if index.I[1] == index.I[2]
            A.value[index] = 0
        end
    end
    return A
end

function enforcesymmetryaroundaxis(u::Field{S})::Field{S} where {S}
    for index in CartesianIndices(u.value)
        if index.I[2] > index.I[1]
            u.value[index] = u.value[index.I[1], index.I[2]]
        end
    end
    return u
end

function J(a::Field{S}, η::Field{S})::NTuple{4, Operator{S}} where {S}

    L1Δη = η*DU*DV + (DV*η)*DU + (DU*η)*DV + (DU*(DV*η))*I
    L1Δa = (a/2)*I
    L2Δη = (1/η)*DU*DV - (1/η^2)*(DU*(DV*η))*I
    L2Δa = (1/a)*DU*DV - (1/a^2)*(DV*a)*DU - (1/a^2)*(DU*a)*DV - (1/a^2)*(DU*(DV*a))*I + (2/a^3)*(DU*a)*(DV*a)*I
    D    = DV + DU
    return ((I-B)*L1Δa + B,  (I-B)*((I-A)*renormalize(L1Δη) + A*D),
            (I-B)*L2Δa,      (I-B)*((I-A)*renormalize(L2Δη) + A*D) + B)
end

function F(a::Field{S}, η::Field{S})::NTuple{2, Field{S}} where {S}
    FonAxis = DV*η + DU*η 
    F1 = η*(DU*(DV*η)) + (DU*η)*(DV*η) + (1/4)*(a^2)
    F2 = (1/a)*(DU*(DV*a)) - (1/a^2)*(DU*a)*(DV*a) + (1/η)*(DU*(DV*η))
    return (enforcesymmetryaroundaxis((I-A)*renormalize(F1) + A*FonAxis),
            enforcesymmetryaroundaxis((I-A)*renormalize(F2) + A*FonAxis))
end

function residual(a::Field{S}, η::Field{S})::NTuple{2, Field{S}} where {S}
    function BCs(F1::Field{S}, F2::Field{S})::NTuple{2, Field{S}} where {S}
        BF1 = (I-B)*F1 + B*(a-bnda) 
        BF2 = (I-B)*F2 + B*(η-bndη) 
        return (BF1, BF2)
    end
    return BCs(F(a, η)...)
end

PS = ProductSpace(ChebyshevGL{U, 16, Float64}(0, 1), 
                  ChebyshevGL{V, 16, Float64}(0, 1))

DU, DV = derivative(PS)
B = incomingboundary(PS)
A = axisboundary(PS)
I = identity(PS)

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


# @test JfromAD ≈ Janalytic
# exit()

using NLsolve
as, ηs = reshapeToTuple(PS, nlsolve(f!, reshapeFromTuple((1 + a, η)); autodiff=:forward, show_trace=true, ftol=1e-9, iterations=100).zero)
# as, ηs = reshapeToTuple(PS, nlsolve(f!, j!, reshapeFromTuple((0.1 + a, η));  show_trace=true, ftol=1e-9, iterations=200).zero)
@show L2(as - a)
@show L2(ηs - η)

using PyPlot
using LaTeXStrings
fig = figure(figsize=(24, 6))
subplot(121) # Create the 1st axis of a 2x2 arrax of axes
PyPlot.title(L"a")
pcolormesh(as)
subplot(122) # Create the 1st axis of a 2x2 arrax of axes
PyPlot.title(L"$\eta$")
pcolormesh(ηs)
show()
