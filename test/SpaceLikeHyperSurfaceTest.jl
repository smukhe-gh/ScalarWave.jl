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
    elseif kind == :initialspacelike
        B[1, :] .= 1
        B[:, 1] .= 1
        return ProductSpaceOperator(space, reshape(diagm(0=>vec(B)), (length(S2), length(S1), length(S2), length(S1))))
    elseif kind == :initialspacelikederivative
        B[2, :] .= 1
        B[:, 2] .= 1  # where do you want to set the derivative conditions? 
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

# u = solve(L, boundary(:initialdata, S)*u0 + boundary(:initialdataderivative, S)*du0)  

# @show maximum(abs(u - u0))
# @show norm(vec(abs(u-u0)))

# FIXME: The problem might very well be with the condition number of the matrix.
# @show cond(L)

# contourf(abs(u-u0), 100)
# using PyPlot
# colorbar()
# show()

#---------------------------------------------------------------
# If I now try using two spacelike hypersurfaces
#---------------------------------------------------------------

struct X1 end
struct X2 end

S = ProductSpace{GaussLobatto{X1, 14,  0,  -1},
                 GaussLobatto{X2, 24,  1,  0}}

dX, dT = derivative(S)
g = Metric{_dd, 2}([1/4, -3/4, 1,4])
D = [dX, dT]
L = sum(g[i,j]*D[i]*D[j] for i in 1:2, j in 1:2)

# replace rows of the operator to impose boundary conditions
L = ⊙(L, boundary(:initialspacelike, S))
L = ⊙(L, boundary(:initialspacelikederivative, S), dT)

# FIXME: The problem might very well be with the condition number of the matrix.
@show cond(L)
@show sort(abs.((eigvals(vec(L)))))

# initial data and solve
u0  = Field(S, (t,x) -> sin(t-x))
du0 = Field(S, (t,x) -> cos(t-x))

u = solve(L, boundary(:initialspacelike, S)*u0 + boundary(:initialspacelikederivative, S)*du0)  

@show maximum(abs(u - u0))
@show norm(vec(abs(u - u0)))

contourf(abs(u-u0), 100)
using PyPlot
colorbar()
show()
