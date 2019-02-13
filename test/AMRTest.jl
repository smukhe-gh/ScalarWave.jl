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

nest1D = conductor(S, condition, 9)
show(nest1D)
exit()

#--------------------------------------------------------------------
# Now try the same with 2D grids
#--------------------------------------------------------------------

# function refinespace(space::Type{S}) where {S<:ProductSpace{GaussLobatto{TagV, NV, maxV, minV}, 
                                                    # GaussLobatto{TagU, NU, maxU, minU}}} where {TagV, NV, maxV, minV, 
                                                                                                # TagU, NU, maxU, minU}

    # SLL = ProductSpace{GaussLobatto{TagV, NV, maxV, (maxV + minV)/2},
                       # GaussLobatto{TagU, NU, maxU, (maxU + minU)/2}}

    # SRR = ProductSpace{GaussLobatto{TagV, NV, (maxV + minV)/2, minV}, 
                       # GaussLobatto{TagU, NU, (maxU + minU)/2, minU}} 

    # SLR = ProductSpace{GaussLobatto{TagV, NV, (maxV + minV)/2, minV},
                       # GaussLobatto{TagU, NU, maxU, (maxU + minU)/2}}

    # SRL = ProductSpace{GaussLobatto{TagV, NV, maxV, (maxV + minV)/2}, 
                       # GaussLobatto{TagU, NU, (maxU + minU)/2, minU}} 

    # return Dict{Array{Int64,1}, Union{Any, Dict}}([1,1]=>SLL, [1,2]=>SLR, 
                                                  # [2,1]=>SRL, [2,2]=>SRR)
# end

# function refinespace(space::Type{S}, threshold::Float64) where {S<:ProductSpace{GaussLobatto{TagV, NV, maxV, minV}, 
                                                                                # GaussLobatto{TagU, NU, maxU, minU}}} where {TagV, NV, maxV, minV, 
                                                                                                                            # TagU, NU, maxU, minU}
    # if prod(maximum(space)) > threshold 
        # return refinespace(S)
    # else
        # return S
    # end
# end

# function prune2D!(nest::Dict, threshold::Float64)
    # for (key, value) in nest
        # if typeof(value) <: Dict
            # prune2D!(value, threshold)
        # else
            # (prod(maximum(value)) > threshold) ? delete!(nest, key) : 0
        # end
    # end
    # return nest
# end

# function driver(S::Type{ProductSpace{GaussLobatto{TagV, NV, maxV, minV},
                                         # GaussLobatto{TagU, NU, maxU, minU}}}, threshold::Float64,
                                                                               # maxlevel::Int) where {TagV, NV, maxV, minV, 
                                                                                                     # TagU, NU, maxU, minU}
    # newdict = Dict{Array{Int64,1}, Any}([1,1]=>S)
    # for level in 1:maxlevel
        # for (key, value) in newdict
            # newdict[key] = refinespace(value, threshold)
        # end
    # end
    # prune2D!(newdict, threshold)
    # return newdict
# end

# nest2D = driver(SUV, 1.0, 6)
# # showdict(nest2D)

# using PyPlot

# function PyPlot. plot(u::Type{S}) where {S<:ProductSpace{GaussLobatto{Tag1, N1, max1, min1},
                                                       # GaussLobatto{Tag2, N2, max2, min2}}} where {Tag1, N1, max1, min1,
                                                                                                   # Tag2, N2, max2, min2}
    # vlines(x=[min2, max2], ymin=min1, ymax=max1)
    # hlines(y=[min1, max1], xmin=min2, xmax=max2)
    # return 0
# end

# function PyPlot. plot(nest::Dict)
    # for (key, value) in nest
        # plot(value)
    # end
    # return 0
# end

# plot(nest2D)
# show()
