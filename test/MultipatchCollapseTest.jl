#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2019
# Simulate Minkowski spacetime on axis
#--------------------------------------------------------------------

using NLsolve, ForwardDiff

function background(PS::ProductSpace{S1,S2})::NTuple{3, Field{ProductSpace{S1, S2}}} where {S1, S2}
    a = Field(PS, (u,v)->1)
    η = Field(PS, (u,v)->(v-u)/2)
    ϕ = Field(PS, (u,v)->1e-2*sin(u+v))
    return (a, η, ϕ)
end

function einstein(boundarydata::NTuple{3, Field{ProductSpace{S1, S2}}}, 
                  initialguess::NTuple{3, Field{ProductSpace{S1, S2}}})::NTuple{3, Field{ProductSpace{S1, S2}}} where {S1, S2}

    function C(a::Field{S}, η::Field{S}, ϕ::Field{S})::NTuple{2, Field{S}} where {S}
        C1 = DU*(DU*η) - (2/a)*(DU*a)*(DU*η) + (4pi*η)*(DU*ϕ)^2
        C2 = DV*(DV*η) - (2/a)*(DV*a)*(DV*η) + (4pi*η)*(DV*ϕ)^2
        return (C1, C2)
    end

    function F(a::Field{S}, η::Field{S}, ϕ::Field{S})::NTuple{3, Field{S}} where {S}
        F1 = DU*(DV*a) - (1/a)*(DU*a)*(DV*a) + (a/η)*(DU*(DV*η)) + (4pi*a)*(DU*ϕ)*(DV*ϕ)
        F2 = DU*(DV*η) + (1/η)*(DU*η)*(DV*η) + (1/4)*(1/η)*(a^2)
        F3 = DU*(DV*ϕ) + (1/η)*(DU*η)*(DV*ϕ) + (1/η)*(DV*η)*(DU*ϕ)

        F1onAxis = DU*DV*a - (1/a)*(DU*a)*(DV*a) + (4pi*a)*(DU*ϕ)*(DV*ϕ)
        F2onAxis = DU*DV*η 
        F3onAxis = DU*DV*ϕ 

        # @show range(a.space.S1)
        # @show range(a.space.S2)
        # @show range(a.space.S1) .== range(a.space.S2)
        # @show all(range(a.space.S1) .== range(a.space.S2))
        if a.space.S1.min == 0 
            return (mix!(mix!(F1, A, F1onAxis), B, a-bnda), 
                    mix!(mix!(F2, A, F2onAxis), B, η-bndη),
                    mix!(mix!(F3, A, F3onAxis), B, ϕ-bndϕ))
        else
            return (mix!(F1, B, a-bnda), 
                    mix!(F2, B, η-bndη),
                    mix!(F3, B, ϕ-bndϕ))
        end
    end

    function f!(fvec::Array{T,1}, x::Array{T,1}) where {T}
        fvec[:] = reshapeFromTuple(F(reshapeToTuple(PS, x)...))
    end

    
    function j!(jvec::Array{T,2}, x::Array{T,1}) where {T}
        jvec[:, :] = ForwardDiff.jacobian(f!, similar(x), x)
    end

    PS = initialguess[1].space
    DU, DV = derivative(PS)
    B = incomingboundary(PS)
    A = axisboundary(PS)

    (bnda, bndη, bndϕ) = combineUVboundary.(computeUboundary(boundarydata), 
                                            computeVboundary(boundarydata), :incoming)

    @show range(PS)

    solved = reshapeToTuple(PS, nlsolve(f!, reshapeFromTuple(initialguess); method=:trust_region, autodiff=:forward,
                                            show_trace=true, ftol=1e-8, iterations=40).zero)
    @show L2.(C(solved...))
    # @show L2.(F(solved...))
    return solved
end

#--------------------------------------------------------
# Compute on many patches
#--------------------------------------------------------

grid = ScalarWave.Grid(Float64, (14, 14), (10,10), (0,1), (0,1))
distribute(grid, background, einstein)
