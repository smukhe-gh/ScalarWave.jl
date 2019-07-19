#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 07-2019
# Test derivatives when the function has a singularity, but 
# not on the grid.
#--------------------------------------------------------------------

struct U end

function compute(N, order)
    S = ChebyshevGL{U, N, Float64}(-1, 1)
    D = derivative(S)
    x = Field(S, x->x)
    return L1(D*(x^order) - order*x^(order-1))
end

# Choose even points to avoid zero. Even then why doesn't
# the derivative converge?
for N in 1:15
    @show N, compute(2N, -1)
end
