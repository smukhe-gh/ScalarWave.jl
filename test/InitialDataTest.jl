#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2019
# Simulate scalar field collpase on axis
#--------------------------------------------------------------------

# [1] Construct Schwarzschild spacetime and evaluate the constraints. 
# [2] Construct rescaled variables and test the rescaled constraint equations
# [3] Construct operators from the rescaled constraints and test using arbitrary functions. 
# [4] Solve for s or ψ, and compare solutions with Schwarzschild.
# [5] Solve for an unknown spacetime. 

function C(a::Field{S}, η::Field{S}, ϕ::Field{S})::NTuple{2, Field{S}} where {S}
    C1 = DU*(DU*η) - (2/a)*(DU*a)*(DU*η) + (4pi*η)*(DU*ϕ)^2
    C2 = DV*(DV*η) - (2/a)*(DV*a)*(DV*η) + (4pi*η)*(DV*ϕ)^2
    return (C1, C2)
end

#===========================================================================#
PS = ProductSpace(ChebyshevGL{U, 18, Float64}(-4, -2), 
                  ChebyshevGL{V, 18, Float64}( 4,  6))
#===========================================================================#
DU, DV = derivative(PS)
I = identity(PS)

# [1] Construct Schwarzschild spacetime and test constraints
M = 1
η = Field(PS, (u,v)->find_r_of_UV(u ,v, M))
f = (32*M^3/η)*exp(-η/2M)
a = sqrt(f)
ϕ = Field(PS, (u,v)->0)
@test all(L2.(C(a, η, ϕ)) .< 1e-10)

# [2] Construct rescaling functions 
function rescale(a::Field{S}, η::Field{S}, ϕ::Field{S})::NTuple{3, Field{S}} where {S}
    r = Field(a.space, (u,v)->(v-u))
    f = a^2
    ψ = sqrt(η/r)
    s = log(f/(2*ψ^4))/(2*r)
    return (s, ψ, ϕ)
end

function unscale(s::Field{S}, ψ::Field{S}, ϕ::Field{S})::NTuple{3, Field{S}} where {S}
    r = Field(s.space, (u,v)->(v-u))
    f = 2*(ψ^4)*exp(2*r*s)
    a = sqrt(f)
    η = r*ψ^2 
    return (a, η, ϕ)
end

# Test if rescaling functions work as they should
@test all(L2.((a, η, ϕ) .- unscale(rescale(a, η, ϕ)...)) .< 1e-12)

# Now check if the constraint equations you found using Mathematica do the right thing. 
function rescaledC(s::Field{S}, ψ::Field{S}, ϕ::Field{S})::NTuple{2, Field{S}} where {S}
    r = Field(s.space, (u,v)->(v-u))
    C1 = (4*π*(DU*ϕ)^2 - 4*(r^2)*ψ*(DU*s)*(DU*ψ) 
        - 2*s*ψ*(DU*r)*(ψ*(DU*r) + 2*r*(DU*ψ)) 
        + (ψ^2)*(DU*(DU*r)) - 2*r*(ψ^2*(DU*r)*(DU*s) + 3*(DU*ψ)^2 - ψ*(DU*(DU*ψ))))
    C2 = (4*π*(DV*ϕ)^2 - 4*(r^2)*ψ*(DV*s)*(DV*ψ) 
        - 2*s*ψ*(DV*r)*(ψ*(DV*r) + 2*r*(DV*ψ)) 
        + (ψ^2)*(DV*(DV*r)) - 2*r*(ψ^2*(DV*r)*(DV*s) + 3*(DV*ψ)^2 - ψ*(DV*(DV*ψ))))
    return (C1, C2)
end

# Test if the r coordinate is being set properly
r = Field(PS, (u,v)->(v-u))
ronUbnd = extractUboundary(r, :incoming)

function rescaledC2onbnd(s::Field{S}, ψ::Field{S}, ϕ::Field{S})::Field{S} where {S} where {Tag}
    r = ronUbnd
    DV = derivative(s.space)
    return (4*π*(DV*ϕ)^2 - 4*(r^2)*ψ*(DV*s)*(DV*ψ) 
        - 2*s*ψ*(DV*r)*(ψ*(DV*r) + 2*r*(DV*ψ)) 
        + (ψ^2)*(DV*(DV*r)) - 2*r*(ψ^2*(DV*r)*(DV*s) + 3*(DV*ψ)^2 - ψ*(DV*(DV*ψ))))
end

@test all(L2.(rescaledC(rescale(a, η, ϕ)...)) .< 1e-10)

# [3] Now make operators out of the constraint functions
function rhs(ψ::Field{S}, ϕ::Field{S})::Field{S} where {S}
    r = Field(ψ.space, (u,v)->(v-u))
    return (4*π*(DU*ϕ)^2 - 6*r*(DU*ψ)^2 + (ψ^2)*(DU*(DU*r)) + 2*r*ψ*(DU*(DU*ψ)))
end

function sOp(ψ::Field{S}, ϕ::Field{S})::Operator{S} where {S}
    r = Field(ψ.space, (u,v)->(v-u))
    L = -2*(ψ^2)*((DU*r)^2)*I - 2*r*(ψ^2)*(DU*r)*DU - 4*r*ψ*(DU*r)*(DU*ψ)*I - 4*(r^2)*ψ*(DU*ψ)*DU
    return L
end

(s, ψ, ϕ) = rescale(a, η, ϕ)
(sglobal, ψglobal, ϕglobal) = (s, ψ, ϕ)
A = sOp(ψ, ϕ)
b = rhs(ψ, ϕ)
@test L2(A*s + b) < 1e-10
@test L2(rescaledC2onbnd(extractUboundary.((s, ψ, ϕ), :incoming)...)) < 1e-10

# Upgrade the functions to line functions
function rhsonbnd(ψ::Field{S}, ϕ::Field{S})::Field{S2} where {S<:ProductSpace{S1, S2}} where {S1, S2}
    r = Field(ψ.space, (u,v)->(v-u))
    r = extractUboundary(r, :incoming)
    ψ = extractUboundary(ψ, :incoming)
    ϕ = extractUboundary(ϕ, :incoming)
    DV = derivative(r.space)
    return (4*π*(DV*ϕ)^2 - 6*r*(DV*ψ)^2 + (ψ^2)*(DV*(DV*r)) + 2*r*ψ*(DV*(DV*ψ)))
end

function sOponbnd(ψ::Field{S}, ϕ::Field{S})::Operator{S2} where {S<:ProductSpace{S1, S2}} where {S1, S2}
    r = Field(ψ.space, (u,v)->(v-u))
    r = extractUboundary(r, :incoming)
    ψ = extractUboundary(ψ, :incoming)
    ϕ = extractUboundary(ϕ, :incoming)
    DV = derivative(r.space)
    I  = identity(r.space)
    L = -2*(ψ^2)*((DV*r)^2)*I - 2*r*(ψ^2)*(DV*r)*DV - 4*r*ψ*(DV*r)*(DV*ψ)*I - 4*(r^2)*ψ*(DV*ψ)*DV
    return L
end

A = sOponbnd(ψ, ϕ)
b = rhsonbnd(ψ, ϕ)
@test L2(A*extractUboundary(s, :incoming) + b) < 1e-10

# Now, compute the nonlinear operator for ψ
function residualψonbnd(ψ::Field{S})::Field{S2} where {S<:ProductSpace{S1, S2}} where {S1, S2}
    r = Field(ψ.space, (u,v)->(v-u))
    r = extractUboundary(r, :incoming)
    s = extractUboundary(sglobal, :incoming)
    ψ = extractUboundary(ψglobal, :incoming)
    ϕ = extractUboundary(ϕglobal, :incoming)
    DV = derivative(r.space)
    I = identity(r.space)
    B = incomingboundary(r.space) + outgoingboundary(r.space)
    ψbnd = B*ψ
    F = (4*π*(DV*ϕ)^2 - 4*(r^2)*ψ*(DV*s)*(DV*ψ) 
       - 2*s*ψ*(DV*r)*(ψ*(DV*r) + 2*r*(DV*ψ)) 
       + (ψ^2)*(DV*(DV*r)) - 2*r*(ψ^2*(DV*r)*(DV*s) + 3*(DV*ψ)^2 - ψ*(DV*(DV*ψ))))
    return (I - B)*F + B*(ψ - ψbnd)
end
 
f = residualψonbnd(ψ)
@test L2(f) < 1e-10

# [4]  Now try solving for the initial data; start with Schwarzschild. 
function sinitialdata(sbnd::Field{S})::Field{S} where {S<:Space{Tag}} where {Tag}
    A = sOponbnd(ψ, ϕ)
    b = rhsonbnd(ψ, ϕ)
    # @show cond(A ⊕ incomingboundary(sbnd.space))
    c = solve(A ⊕ incomingboundary(sbnd.space), sbnd ⊕ (-b))
    return c
end

sonUbnd = extractUboundary(s, :incoming)
sbnd = incomingboundary(sonUbnd.space)*sonUbnd

ssolved = sinitialdata(sbnd)
@test L2(ssolved - sonUbnd) < 1e-10
@test L2(rescaledC2onbnd(ssolved, extractUboundary(ψ, :incoming),
                         extractUboundary(ϕ, :incoming))) < 1e-10

# Compute ψ for Schwarzschild with the nonlinear solver 
using NLsolve, ForwardDiff
ψonUbnd = extractUboundary(ψ, :incoming)

function Fforψ(ψ::Field{S})::Field{S} where {S}
    r = ronUbnd
    s = extractUboundary(sglobal, :incoming)
    ϕ = extractUboundary(ϕglobal, :incoming)
    DV = derivative(s.space)
    I  = identity(s.space)
    B  = incomingboundary(s.space) + outgoingboundary(s.space)
    ψbnd = B*ψonUbnd

    F = (4*π*(DV*ϕ)^2 - 4*(r^2)*ψ*(DV*s)*(DV*ψ) 
       - 2*s*ψ*(DV*r)*(ψ*(DV*r) + 2*r*(DV*ψ)) 
       + (ψ^2)*(DV*(DV*r)) - 2*r*(ψ^2*(DV*r)*(DV*s) + 3*(DV*ψ)^2 - ψ*(DV*(DV*ψ))))
    return (I - B)*F + B*(ψ - ψbnd)
end

function f!(fvec::Array{T,1}, x::Array{T,1}) where {T}
    fvec[:] = reshape(Fforψ(reshape(sonUbnd.space, x)))
end

function j!(jvec::Array{T,2}, x::Array{T,1}) where {T}
    jvec[:, :] = ForwardDiff.jacobian(f!, similar(x), x)
end

# ψsolved = reshapeToTuple(PS, nlsolve(f!, j!,  reshapeFromTuple((a0, η0, ϕ0)); method=:trust_region, show_trace=true, ftol=1e-10, iterations=120).zero)
ψsolved = reshape(ψonUbnd.space, nlsolve(f!, reshape(ψonUbnd/ψonUbnd); method=:trust_region, autodiff=:forward, show_trace=true, ftol=1e-10, iterations=120).zero)
@test L2(ψsolved - ψonUbnd) < 1e-10

# [5] Now solve with arbitrary initial data and understand how they behave.
# First solve for s
snew = Field(PS, (u,v)->(v-u))
ψnew = Field(PS, (u,v)->1)
p = 0.1
ϕnew = Field(PS, (u,v)->p*(exp(-(u-1/2)^2) + exp(-(v-1/2)^2)))

function sinitialdatanew(sbnd::Field{S})::Field{S} where {S<:Space{Tag}} where {Tag}
    A = sOponbnd(ψnew, ϕnew)
    b = rhsonbnd(ψnew, ϕnew)
    # @show cond(A ⊕ incomingboundary(sbnd.space))
    c = solve(A ⊕ incomingboundary(sbnd.space), sbnd ⊕ (-b))
    return c
end

sonUbnd = extractUboundary(snew, :incoming)
sbnd = incomingboundary(sonUbnd.space)*sonUbnd

ssolved = sinitialdatanew(sbnd)
@test L2(ssolved - sonUbnd) < 1e-10
@test L2(rescaledC2onbnd(ssolved, extractUboundary(ψnew, :incoming),
                         extractUboundary(ϕnew, :incoming))) < 1e-10
