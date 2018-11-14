#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
# Compute error estimates
#--------------------------------------------------------------------
# TODO: Compute L2 norm
# TODO: Compute relative L2 norm error
# TODO: Compute Linf norm
# TODO: Compute error estimate from coefficents


function L2normError(A::Field{S}, B::Field{S})::Real
    D = (A - B)^2

end
