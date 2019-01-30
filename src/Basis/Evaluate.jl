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

function (c::Field{S})(x::Number)::Number where {S<:Chebyshev{Tag, N, max, min}} where {Tag, N, max, min}
    @assert x in c.space
    xmap = (x - (max + min)/2)*(2/(min - max))
    return - sum(prefactor(order+1, N)*cheb(order, xmap)*c.value[order+1] for order in 0:N)
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
