#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2019
# Test basis transformations
#--------------------------------------------------------------------

S = ChebyshevGL{U, 6, Float64}(-1, 1)
f = Field(S, x->x)

@show basistransform(f)
@show basistransform(f^2)
@show basistransform(f^3)
@show basistransform(f^4)

