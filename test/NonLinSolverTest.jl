#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2019
# Test a nonlinear solver using Newton iterations
#--------------------------------------------------------------------

# Solve a simple non-linear BVP problem for confidence
# u_xx = Exp[u]; u(+/-1) =  0

# set up space
struct X end
S = GaussLobatto(X, 16)
x = Field(S, x->x)
D = derivative(S)
I = eye(S)
B = boundary(S)

# start with an initial guess
uold = Field(S, x->0) 

tol = 1e-14
maxiter = 300

# now iterate [setting up the iteration can be non-trivial]
function iterate(tol, maxiter::Int)
    for iteration in range(1, length=maxiter)
        unew  = solve(D*D + B, exp(uold)) 
        error = maximum(abs(unew - uold)) 
        @show iteration, error
        if error < tol 
            return unew
        else
            for index in eachindex(uold.value)
                uold.value[index] = unew.value[index]
            end
        end
    end
end

unew = iterate(tol, maxiter)
println("u(0) = ", unew.value[Int(order(S)/2)+1])
println("===============================================================")
println("===============================================================")
println("===============================================================")
println("===============================================================")

# The above problem doesn't give much intuition as to where the iterations
# are coming from and why they would work. 
# Below we set up an extensible model for solving general non-linear
# partial differential equation: L(u) = F(u) where L is the principle part. 
# We define two operators F and Fu
# which generates a linear iteration for delta, i.e. u(i+1) = u(i) + delta

function F(u::Field)
    return exp(u)
end

function Fu(u::Field)
    return exp(u)
end

# start with an initial guess
uold = Field(S, x->0) 

function iterateNewtonKantorovich(F::Function, Fu::Function, tol::Float64, maxiter::Int)
    for iteration in range(1, length=maxiter)
        if iteration == maxiter
            @warn "Iterations did not converge."
            return 0
        end
        delta = solve(D*D - Fu(uold)*I,  D*D*uold - F(uold)) 
        unew  = uold + delta 
        error = maximum(abs(unew - uold)) 
        @show iteration, error
        if error < tol 
            return unew
        else
            for index in eachindex(uold.value)
                uold.value[index] = unew.value[index]
            end
        end
    end
end

unew = iterateNewtonKantorovich(F, Fu, tol, maxiter)
