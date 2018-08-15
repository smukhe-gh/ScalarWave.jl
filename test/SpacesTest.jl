#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2018
#--------------------------------------------------------------------

"""
We want to be able to define two operations at 
the moment.
    -- Multiplication of two fields
    -- Action of a derivative operator on a field in 1D.
We assume all our datastructures are ordered sets
"""

# Define a manifold
M = Manifold(-1, 1, 10)

# Define a chart over the manifold
S = chart(M, "Chebyshev")

# Now, given a chart, compute the field 
# over manifold. Note that the field doesn't
# need to know about the underlying chart. 
u = Field(M, map(x->x^2+3, S)) 

# Define another field on the manifold and multiply the
# two fields together
v = Field(M, map(x->x^3, S))
w = u*v

# For taking the derivative, promote the 1D field to a vector
# NOTE: Ask Erik: is this reasonable?
uvec  = VectorSpace(u, M.npoints, u.elements)

# Let's define it's dual
udual = DualSpace(uvec, u.elements.^-1) 

# Take the inner product to get the dimension of the vectorspace
# NOTE: Do we need a norm?
scalar = udual*uvec

# Now we define the derivative operator
D = derivOperator(u)

# Compute the action of the derivative operator on u
u' = D*u
