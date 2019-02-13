#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 02-2019
#--------------------------------------------------------------------

function refine(u::Field{S}) where {S<:GaussLobatto{Tag, N, max, min}} where {Tag, N, max, min}
    SL = GaussLobatto{Tag, N, max, (max + min)/2} 
    SR = GaussLobatto{Tag, N, (max + min)/2, min} 
    uL = Field(SL, x->u(x))
    uR = Field(SR, x->u(x))
    return Dict{Array{Int64,1}, Union{Field, Dict}}([1]=>uL, [2]=>uR)
end

function refine(u::Field{S}) where {S<:ProductSpace{GaussLobatto{TagV, NV, maxV, minV}, 
                                                    GaussLobatto{TagU, NU, maxU, minU}}} where {TagV, NV, maxV, minV, 
                                                                                                TagU, NU, maxU, minU}

    SLL = ProductSpace{GaussLobatto{TagV, NV, maxV, (maxV + minV)/2},
                       GaussLobatto{TagU, NU, maxU, (maxU + minU)/2}}

    SRR = ProductSpace{GaussLobatto{TagV, NV, (maxV + minV)/2, minV}, 
                       GaussLobatto{TagU, NU, (maxU + minU)/2, minU}} 

    SLR = ProductSpace{GaussLobatto{TagV, NV, (maxV + minV)/2, minV},
                       GaussLobatto{TagU, NU, maxU, (maxU + minU)/2}}

    SRL = ProductSpace{GaussLobatto{TagV, NV, maxV, (maxV + minV)/2}, 
                       GaussLobatto{TagU, NU, (maxU + minU)/2, minU}} 

    uLL = Field(SLL, (x,y)->u(x,y))
    uRR = Field(SRR, (x,y)->u(x,y))
    uRL = Field(SRL, (x,y)->u(x,y))
    uLR = Field(SLR, (x,y)->u(x,y))

    return Dict{Array{Int64,1}, Union{Field, Dict}}([1,1]=>uLL, [1,2]=>uLR, 
                                                    [2,1]=>uRL, [2,2]=>uRR)
end

function refine(space::Type{S}) where {S<:GaussLobatto{Tag, N, max, min}} where {Tag, N, max, min}
    SL = GaussLobatto{Tag, N, max, (max + min)/2} 
    SR = GaussLobatto{Tag, N, (max + min)/2, min} 
    return Dict{Array{Int64,1}, Union{Any, Dict}}([1]=>SL, [2]=>SR)
end

# function refine(u::Type{S}) where {S<:ProductSpace{GaussLobatto{TagV, NV, maxV, minV}, 
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

    # return Dict{Array{Int64,1}, Union{S, Dict}}([1,1]=>SLL, [1,2]=>SLR, 
                                                # [2,1]=>SRL, [2,2]=>SRR)
# end

