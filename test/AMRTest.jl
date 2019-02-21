#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2019
# Plot fields using PyPlot
#--------------------------------------------------------------------

S   = GaussLobatto{U, 100, 3.0, -3.0}
SUV = ProductSpace{GaussLobatto{V, 20,  5.0, -5.0},
                   GaussLobatto{U, 20,  5.0, -5.0}}

y = Field(S, x->sin(x))
ϕ = Field(SUV, (U,V)->sinpi(U)*cospi(V))
yLRdict = refine(y)
ϕLRdict = refine(ϕ) 

@testset "AMR" begin
    @test length(yLRdict) == 2
    @test length(ϕLRdict) == 4
end

function condition(S::Type{GaussLobatto{Tag, N, max, min}})::Bool where {Tag, N, max, min}
    maximum(S) > 1 ? (return 1) : (return 0) 
end

# nest1D = conductor(S, condition, 9)
# show(nest1D)
# exit()

function condition2D(S::Type{ProductSpace{S1, S2}})::Bool where {S1, S2}
    prod(maximum(S)) > 1 ? (return 1) : (return 0) 
end

function condition2Dexcise(S::Type{ProductSpace{S1, S2}})::Bool where {S1, S2}
    UV = Field(ProductSpace{S1, S2}, (x,y)->(x-y))
    for value in vec(UV)
        if abs(value) < 0.2
            return 1
        end
    end
    return 0
end

# nest2D = conductor(SUV, condition2Dexcise, 5)
# plot(nest2D)
# using PyPlot
# axis("equal")
# show()

function Field(dictionary::Dict{Array{Int64,1}, Union{Any, Dict}}, map::Function)
    for (key, value) in dictionary
        dictionary[key] = Field(value, map)
    end
    return dictionary
end
    
nest2DField = Field(nest2D, (x,y)->(x+y))
levels(nest2DField, globallength=20)
# contourf(nest2DField, 100, globalmax=maximum(nest2DField), globalmin=minimum(nest2DField), globallevels=levels(nest2DField, globallength=20))
# show()
    
