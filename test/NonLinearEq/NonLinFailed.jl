#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 07-2019
# Test non-linear solver
#--------------------------------------------------------------------

using NLsolve

struct U end
S = ChebyshevGL{U, 22, Float64}(-1,1)
D = derivative(S)
I = identity(S)
B = incomingboundary(S) ⊕ outgoingboundary(S)

function Base. exp(u::Field{S}) where {S}
    return Field(u.space, exp.(u.value))
end

function f!(F, x)
    u = Field(S, x)
    F[:] = reshape((I-B)*(D*D*u - exp(u)))
end

function j!(J, x)
    u = Field(S, x)
    J[:, :] = reshape(D*D - exp(u)*I)
end

u0 = reshape(Field(S, u->(u-1)*(u+1)))
u  = nlsolve(f!, j!, u0, method=:newton, show_trace=true, ftol=1e-10)
@show converged(u)

using PyPlot
x = Field(S, x->x)
plot(x.value, u.zero)
show()
