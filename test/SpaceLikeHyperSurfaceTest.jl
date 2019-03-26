#---------------------------------------------------------------
# Test evolution using a (pair of) spacelike hypersurface (s)
# Soham M 3/2019
#---------------------------------------------------------------

using LinearAlgebra

struct X end
struct T end

S = ProductSpace{GaussLobatto{X, 14, pi, -pi},
                 GaussLobatto{T, 24,  1,  0}}
dX, dT = derivative(S)

# Solve the wave equation in 1D, with Newmann and Dirichlet
# boundary conditions at the surface
L = dT*dT - dX*dX

# construct your boundary operator
function boundary(kind::Symbol, space::Type{ProductSpace{S1, S2}}) where {S1, S2}
    B = zeros(size(space))
    if kind == :initialdata 
        B[1, :] .= 1
        return ProductSpaceOperator(space, reshape(diagm(0=>vec(B)), (length(S2), length(S1), length(S2), length(S1))))
    elseif kind == :initialdataderivative
        B[2, :] .= 1  # where do you want to set the derivative conditions? 
        return ProductSpaceOperator(space, reshape(diagm(0=>vec(B)), (length(S2), length(S1), length(S2), length(S1))))
    else
        @warn "Invalid kind. Returning identity operator"
        return eye(space)
    end
end

# replace rows of the operator to impose boundary conditions
L = ⊙(L, boundary(:initialdata, S))
L = ⊙(L, boundary(:initialdataderivative, S), dT)

# initial data and solve
u0  = Field(S, (t,x) -> sin(t-x))
du0 = Field(S, (t,x) -> cos(t-x))

u = solve(L, boundary(:initialdata, S)*u0 + boundary(:initialdataderivative, S)*du0)  

@show maximum(abs(u - u0))
@show norm(vec(abs(u-u0)))

# FIXME: The problem might very well be with the condition number of the matrix.
@show cond(L)

contourf(abs(u-u0), 100)
using PyPlot
colorbar()
show()
