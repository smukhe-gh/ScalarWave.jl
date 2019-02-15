#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2019
# Test a nonlinear solver using Newton iterations
#--------------------------------------------------------------------

using LinearAlgebra
# u_xx = Exp[u]; u(+/-1) =  0

struct X end
S = GaussLobatto(X, 36)
x = Field(S, x->x)
D = derivative(S)
B = boundary(S)
SI = eye(S)

abstol  = 1e-10
maxiter = 40

function LinearAlgebra. norm(u::Field{S})
    return norm(u.value)
end

function nonlinres(u)
    return D*D*u - exp(u)
end

function linOP(u)
    return D*D - exp(u)*SI
end

function linRHS(u)
    return (D*D)*u - exp(u)
end

u0 = Field(S, x->0) 
b  = Field(S, x->0) 

function iterate(maxiter::Int, abstol::Float64, uguess::Field{S}) where {S}
    u = uguess
    println("Starting non-linear solve")
    for iteration in 1:maxiter
        Δu = solve(linOP(u) + B, (-1)*linRHS(u))
        u  = u + Δu
        error = norm(nonlinres(u))
        println("iter = $iteration, error = $error") 
        (error < abstol) ? (return (u, error)) : 0
    end
    return (u, error)
end

(uf, err) = iterate(maxiter, abstol, u0)
@test err < abstol
