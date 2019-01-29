#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2019
# Evaluate a field function at an arbitrary point it's space
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Start with the 1D case 
#--------------------------------------------------------------------

function Base. in(coordinate::Number, S)::Bool 
    if maximum(S) >= coordinate >= minimum(S)
        return true
    else
        return false
    end
end

function (c::Field{S})(x::Number)::Number where {S<:Chebyshev{Tag, N}} where {Tag, N}
    @assert x in c.space
    return sum(c.value[m]*chebmod(m-1, Float64(x)) for m in 1:(order(S)+1))
end

function (u::Field{S})(x::Number)::Number where {S<:GaussLobatto{Tag, N}} where {Tag, N}
    @assert x in u.space
    c = basistransform(u)
    return c(x)
end

#--------------------------------------------------------------------
# Now, work with the 2D case 
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# The 1D grid case 
#--------------------------------------------------------------------
#
#--------------------------------------------------------------------
# The 2D grid case 
#--------------------------------------------------------------------
