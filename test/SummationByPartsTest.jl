#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2019
# Test summation-by-parts 
#--------------------------------------------------------------------

N = 2

# Test with finite difference operators
S = GaussLobatto{U, N, 1, -1}
D = derivative(S)
E = Operator(S, [((i == j == 1) ? -1 : (i == j == N + 1 ? 1 : 0)) for i in 1:N+1, j in 1:N+1])
W = Operator(S, [((i == j) ? chebw(i,N) : 0) for i in 1:N+1, j in 1:N+1])


@show chebw(1, 2)
@show chebw(2, 2)
@show chebw(3, 2)
@show W 
@show D
@show (W*D).value + transpose((W*D).value)


u = Field(S, x->x^2)
x = Field(S, x->x)
@show  
@show sum((W*u).value)
