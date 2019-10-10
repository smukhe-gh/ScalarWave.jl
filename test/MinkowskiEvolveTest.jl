#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2019
# Simulate Minkowski spacetime on axis
#--------------------------------------------------------------------

using NLsolve

function minkowskisolver(boundarydata::NTuple{2, Field{ProductSpace{S1, S2}}}, 
                         initialguess::NTuple{2, Field{ProductSpace{S1, S2}}})::NTuple{2, Field{ProductSpace{S1, S2}}} where {S1, S2}


    function F(a::Field{S}, ω::Field{S})::NTuple{2, Field{S}} where {S}

        function regularize!(F::Field{S})::Field{S} where {S}
            FonAxis = DV*r + DU*r 
            for index in CartesianIndices(F1)
                if index.I[1] == index.I[2] && r.value[index] <= eps(Float64)
                    F.value[index] = FonAxis.value[index]
                end
            end
            return F
        end

        r  = η*ω
        F1 = r*(DU*(DV*r)) + (DU*r)*(DV*r) + (1/4)*(a^2)
        F2 = (1/a)*(DU*(DV*a)) - (1/a^2)*(DU*a)*(DV*a) + (1/r)*(DU*(DV*r))

        return (regularize!(F1), regularize!(F2))
    end

    function C(a::Field{S}, r::Field{S})::NTuple{2, Field{S}} where {S}
        DU, DV = derivative(a.space)
        C1 = DU*(DU*r) - (2/a)*(DU*a)*(DU*r) 
        C2 = DV*(DV*r) - (2/a)*(DV*a)*(DV*r) 
        return (C1, C2)
    end

    function residual(a::Field{S}, r::Field{S})::NTuple{2, Field{S}} where {S}
        function BCs(F1::Field{S}, F2::Field{S})::NTuple{2, Field{S}} where {S}
            BF1 = (I-B)*F1 + B*(a-bnda) 
            BF2 = (I-B)*F2 + B*(r-bndr) 
            return (BF1, BF2)
        end
        return BCs(F(a, r)...)
    end


    function f!(fvec::Array{T,1}, x::Array{T,1}) where {T}
        fvec[:] = reshapeFromTuple(residual(reshapeToTuple(PS, x)...))
    end

    function reshapeFromTuple(U::NTuple{2, Field})
        return vcat(reshape(U[1]), reshape(U[2]))
    end
    
    function reshapeToTuple(space::S, x::Array{T,1})::NTuple{2, Field}  where {S, T}
        U = reshape(x, :, 2)
        return (reshape(space, U[:, 1]), reshape(space, U[:, 2]))
    end

    η = Field(PS, (u,v)->(v-u)/2)
    (bnda, bndr) = boundarydata
    (a0, r0) = initialguess

    @show L2.(C(initialguess...))
    @show L2.(F(initialguess...))

    solved = reshapeToTuple(PS, nlsolve(f!, reshapeFromTuple((a0, r0)); 
                                            show_trace=true, ftol=1e-9, iterations=50).zero)
    @show L1.(F(solved...))
    @show L2.(C(solved...))
    @show L2.(F(solved...))
    return solved
end

PS = ProductSpace(ChebyshevGL{U, 6, Float64}(0, 1), 
                  ChebyshevGL{V, 6, Float64}(0, 1))

DU, DV = derivative(PS)
B = incomingboundary(PS)
I = identity(PS)

a0 = Field(PS, (u,v)->1)
ω0 = Field(PS, (u,v)->1)

(asol, rsol) = minkowskisolver((B*a0, B*ω0), (1 + a0, ω0))
pcolormesh(asol)
show()
pcolormesh(rsol)
show()

