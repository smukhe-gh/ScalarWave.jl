#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
# Compute error estimates
#--------------------------------------------------------------------
# TODO: Compute L2 norm
# TODO: Compute relative L2 norm error
# TODO: Compute Linf norm
# TODO: Compute error estimate from coefficents

function L2Error(A::Field{S}, B::Field{S})::Real where {S}
    WUV = integral(S)
    L2  = sqrt(abs(WUV*(A - B)^2))
    return L2
end

function L2ErrorRelative(A::Field{S}, B::Field{S})::Real where {S}
    WUV = integral(S)
    L2  = sqrt(abs(WUV*(A - B)^2))/sqrt(abs(WUV*B))
    return L2
end
