#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 07-2019
# Test non-linear solver
#--------------------------------------------------------------------

struct U end
struct V end

SU  = ChebyshevGL{U, 10, Float64}(-1, 1)
SV  = ChebyshevGL{V, 10, Float64}(-1, 1)
SUV = ProductSpace(SU, SV)
DU, DV = derivative(SUV)
B   = incomingboundary(SUV)
I   = identity(SUV)

function stack(u::Field{S}) where {S}
    return reshape(u) 
end

function unstack(space::S, x::Array{T,1}) where {S, T}
    return reshape(space, x)
end

function f!(F, x)
    ϕ = unstack(SUV, x)
    res = (I-B)*(DU*DV*ϕ)  
    F[:] = stack(res)
end

ϕ = Field(SUV, (u,v)->exp(-u^2/0.1) + exp(-v^2/0.1))
noise = Field(SUV, (u,v)->rand())
x0 = stack(ϕ ⊕ noise)

using NLsolve
u = nlsolve(f!, x0, method = :newton, show_trace = true)
@show converged(u)

F = stack(ϕ - ϕ)
f!(F, u.zero)
@show F

# pcolormesh(reshape(SUV, u.zero))
# show()

