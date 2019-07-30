#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 07-2019
# Let's write our own awesome non-linear solver.
# NOTE: We're imposing homogenous boundary conditions on Δ
#--------------------------------------------------------------------

function Base. exp(u::Field{S})::Field{S} where {S}
    return Field(u.space, exp.(u.value))
end

function J(u::Field{S})::Operator{S} where {S}
    L = D*D - exp(u)*I
    return (I-B)*L + B
end

function F(u::Field{S})::Field{S} where {S}
    f = D*D*u - exp(u)
    return (I-B)*f + B*(u-b)
end

function NewtonStep(u::Field{S}, J::Operator{S}, F::Field{S})::Field{S} where {S}
    Δ = solve(J, -F)
    return u + Δ
end

function NewtonIterator(u::Field{S}, J::Function, F::Function; maxiter=1000, abstol=1e-14, reltol=1e-14) where {S}
    println("Initial Newton residual = ", L2(F(u)))
    for iteration in 1:maxiter
        u = NewtonStep(u, J(u), F(u))
        println(iteration, "\t", L2(F(u)))
        L2(F(u)) < abstol ? (return u) : continue
    end
    return u 
end 
    
struct M end
S = ChebyshevGL{M, 15, BigFloat}(-1, 1)
D = derivative(S)
I = identity(S)
B = incomingboundary(S) ⊕ outgoingboundary(S)

s = Field(S, x->1)
b = Field(S, x->1)

u0 = (I-B)*s + B*b 
uf = NewtonIterator(u0, J, F, abstol=1e-20)

# using PyPlot
# @show uf.value[8]
# plot(uf)
# show()
