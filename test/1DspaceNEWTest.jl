#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 04-2019
# Test new functions for nD spaces 
#--------------------------------------------------------------------

function collocation(::Type{S}, i::Int) where {S<:GaussLobatto{Tag, N, max, min}} where {Tag, N, max, min} 
    @assert max > min
    @assert i <= N+1
    return chebx(i, N)*(min - max)/2 + (max + min)/2
end

function collocation(S::Type{T}, index::CartesianIndex) where {T<:GaussLobatto{Tag}} where {Tag}
    return collocation(S, index[1])
end

function collocation(S::Type{ProductSpace{S1, S2}}, index::CartesianIndex)::NTuple{2} where {S1, S2}
    return (collocation(S1, index[1]), collocation(S2, index[2]))
end

function project(S::Type{T}, map::Function)::Field{T} where {T}
    value = zeros(size(S))
    for index in CartesianIndices(value)
        value[index] = map(collocation(S, index)...)
    end
    return Field(S, value)
end

S = GaussLobatto{U, 3, 1, -1}
x = project(S, x->x)
y = project(S, x->1)

@test maximum(S) ==  1
@test minimum(S) == -1
@test order(S)   ==  3
@test size(S)    ==  4
@test identity(S)*x ≈ x
@test derivative(S)*x ≈ y 
@test sum((integral(S)*x).value) < 1e-15
println("Finished 1D tests")

S2 = ProductSpace{GaussLobatto{U, 3, 3, -1},
                  GaussLobatto{V, 5, 1, -1}}

u  = project(S2, (x,y)->x+y)
v  = project(S2, (x,y)->x)
w  = project(S2, (x,y)->y)

@test maximum(S2) == ( 3, 1)
@test minimum(S2) == (-1,-1)
@test order(S2)   == ( 3, 5)
@test size(S2)    == ( 4, 6)


@test identity(S2)*u ≈ u
@test sum((integral(S2)*u).value) ≈ 8.0
@test derivative(S2)[2]*u ≈ w 
@test derivative(S2)[1]*u ≈ v 
println("Finished 2D tests")
