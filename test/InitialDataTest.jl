#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 07-2019
# Test functions inside functions
# TODO: Try with Neumann boundary condition at the end
#--------------------------------------------------------------------

using NLsolve, PyPlot

function solveODE(a::Field{S}, r::Field{S}, ϕ::Field{S})::Field{S} where {S<:Space{Tag}} where {Tag}
    D = derivative(ϕ.space)
    I = identity(ϕ.space)
    B = incomingboundary(ϕ.space) + outgoingboundary(ϕ.space)
    A = D*D - (2/a)*(D*a)*D + 4pi*(D*ϕ)^2*I
    temp = solve(A ⊕ B, B*r)
    @show L2(A*temp)
    return temp
end

function solveODE2(a::Field{S}, r::Field{S}, ϕ::Field{S})::Field{S} where {S<:Space{Tag}} where {Tag}
    D = derivative(ϕ.space)
    I = identity(ϕ.space)
    B = incomingboundary(ϕ.space) ⊕ outgoingboundary(ϕ.space)*D
    A = D*D - (2/a)*(D*a)*D + 4pi*(D*ϕ)^2*I
    L = ⊕(A, incomingboundary(ϕ.space) + outgoingboundary(ϕ.space), B)
    temp = solve(A ⊕ B, B*r)
    @show L2(A*temp)
    return temp
end

function initialdatasolver(PS::ProductSpace{S1, S2}, p::Number)::Field{ProductSpace{S1, S2}} where {S1, S2}
    f0 = Field(PS, (u,v)->1)
    ϕ0 = Field(PS, (u,v)->p*exp(-(v-4)^2) + p*exp(-(u+7)^2))
    r0 = Field(PS, (u,v)->v-u)
    a0 = -sqrt(2*f0) 

    varU = solveODE2(extractUboundary(a0), extractUboundary(r0), extractUboundary(ϕ0))
    varV = solveODE2(extractVboundary(a0), extractVboundary(r0), extractVboundary(ϕ0))
    varC = combineUVboundary(varU, varV)
    return varC
end

#--------------------------------------------------------------------
# Test initial data families
#--------------------------------------------------------------------

struct U end
struct V end

PS = ProductSpace(ChebyshevGL{U, 4, Float64}(-8, -6), 
                  ChebyshevGL{V, 4, Float64}( 3,  5))

# using PyPlot
plot(extractUboundary(Field(PS, (u,v)->v-u)), plotstyle="r--")
plot(extractUboundary(initialdatasolver(PS, 0.1)), plotstyle="-")
show()

