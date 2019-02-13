#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 02-2019
# Refinement routines and drivers for Fields and Spaces
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Refine Fields/Spaces
#--------------------------------------------------------------------

function refine(u::Field{S}) where {S<:GaussLobatto{Tag, N, max, min}} where {Tag, N, max, min}
    SL = GaussLobatto{Tag, N, max, (max + min)/2} 
    SR = GaussLobatto{Tag, N, (max + min)/2, min} 
    uL = Field(SL, x->u(x))
    uR = Field(SR, x->u(x))
    return Dict{Array{Int64,1}, Union{Field, Dict}}([1]=>uL, [2]=>uR)
end

function refine(space::Type{S}) where {S<:GaussLobatto{Tag, N, max, min}} where {Tag, N, max, min}
    SL = GaussLobatto{Tag, N, max, (max + min)/2} 
    SR = GaussLobatto{Tag, N, (max + min)/2, min} 
    return Dict{Array{Int64,1}, Union{Any, Dict}}([1]=>SL, [2]=>SR)
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

function refine(u::Type{S}) where {S<:ProductSpace{GaussLobatto{TagV, NV, maxV, minV}, 
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

    return Dict{Array{Int64,1}, Union{Any, Dict}}([1,1]=>SLL, [1,2]=>SLR, 
                                                [2,1]=>SRL, [2,2]=>SRR)
end

#--------------------------------------------------------------------
# Driver routines 
#--------------------------------------------------------------------

function driver(S::Type{GaussLobatto{Tag, N, max, min}}, condition::Function) where {Tag, N, max, min}
    condition(S) ? (return refine(S)) : (return S)
end

function driver(S::Type{ProductSpace{GaussLobatto{Tag1, N1, max1, min1},
                                     GaussLobatto{Tag2, N2, max2, min2}}}, condition::Function) where {Tag1, N1, max1, min1,
                                                                                                       Tag2, N2, max2, min2}
    condition(S) ? (return refine(S)) : (return S)
end

function driver(dictionary::Dict{Array{Int64,1}, Union{Any, Dict}}, condition::Function)
    for (key, value) in dictionary
        dictionary[key] = driver(value, condition)
    end
    return dictionary
end

function prune!(dictionary::Dict, condition::Function)
    for (key, value) in dictionary 
        if typeof(value) <: Dict
            prune!(value, condition)
        else
            condition(value) ? delete!(dictionary, key) : 0
        end
    end
    return dictionary
end

function Base. show(dictionary::Dict{Array{Int64,1}, Union{Any, Dict}})
    for (key, value) in dictionary
        if typeof(value) <: Dict
            show(value)
        else
            println(value) 
        end
    end
end

function conductor(S::Type{GaussLobatto{Tag, N, max, min}}, condition::Function, maxlevel::Int) where {Tag, N, max, min}
    dictionary = Dict{Array{Int64,1}, Any}([1]=>S)
    for level in 1:maxlevel
        for (key, value) in dictionary
            dictionary[key] = driver(value, condition)
        end
    end
    return prune!(dictionary, condition)
end

function conductor(S::Type{ProductSpace{GaussLobatto{Tag1, N1, max1, min1},
                                        GaussLobatto{Tag2, N2, max2, min2}}}, condition::Function, maxlevel::Int) where {Tag1, N1, max1, min1,
                                                                                                                         Tag2, N2, max2, min2}
    dictionary = Dict{Array{Int64,1}, Union{Any, Dict}}([1,1]=>S)
    for level in 1:maxlevel
        for (key, value) in dictionary
            dictionary[key] = driver(value, condition)
        end
    end
    return prune!(dictionary, condition)
end


