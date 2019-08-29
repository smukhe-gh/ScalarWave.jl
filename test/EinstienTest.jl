#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 07-2019
# Write an initial data solver
# We have the freedom to choose two out of three initial fields
# at the incoming boundaries. This leads to a second order
# differential equation in r, which requires two boundary conditions 
# We specify the starting and final values of r which in turn decides 
# the span of r, and set f = const as a gauge choice.  
# TODO: Solve a more accurate version of the initial data since it's
# cheap.
#--------------------------------------------------------------------

using DualNumbers
using NLsolve
using Roots

function IDSolver(r::Field{S}, ϕ::Field{S})::Field{S} where {S<:Space{Tag}} where {Tag}

    #  Construct local operators
    D = derivative(ϕ.space)
    I = identity(ϕ.space)
    B = incomingboundary(ϕ.space) + outgoingboundary(ϕ.space)

    # Construct RHS to enforce boundary conditions
    b  = B*r 

    # Construct operator after choosing the gauge f = 1
    A = -2*D*D - ((D*ϕ)^2)*I

    # Solve  the 1D equation
    rsolved = solve(A ⊕ B, b)

    return rsolved
end

function extractUboundary(u::Field{ProductSpace{S1,S2}})::Field{S1} where {S1, S2}
    return Field(u.space.S1, u.value[1,:])
end

function extractVboundary(u::Field{ProductSpace{S1,S2}})::Field{S2} where {S1, S2}
    return Field(u.space.S2, u.value[:,1])
end

function combineUVboundary(u::Field{S1}, v::Field{S2})::Field{ProductSpace{S1, S2}} where {S1, S2}
    PS = ProductSpace(u.space, v.space)
    w  = Field(PS, (u,v)->1) 
    w.value[1,:] = u.value
    w.value[:,1] = v.value
    return w
end

#--------------------------------------------------------------------
# Core functions 
#--------------------------------------------------------------------

function f!(fvec::Array{T,1}, x::Array{T,1}) where {T}
    fvec[:] = reshapeFromTuple(F(reshapeToTuple(S, x)...))
end

function j!(J::Array{T,2}, x::Array{T,1}) where {T}
    x = Array{Union{eltype(x), Dual{eltype(x)}}}(x) 
    for index in 1:length(x)
        J[:, index] = Δf(index, x)
    end
    x = Array{eltype(x)}(x) 
end

function Δf(j::Int, x::Array{T,1})::Array{T,1} where {N, T <: Union{X, Dual{X}}} where {X}
    x[j] = Dual(x[j], 1)
    δf   = reshapeFromTuple(F(reshapeToTuple(S, x)...)) 
    x[j] = realpart(x[j])
    return dualpart.(δf)
end

function Base. +(a::Number, u::Field{S})::Field{S} where {S}
    return Field(u.space, a .+ u.value)
end

function Base. -(u::Field{S})::Field{S} where {S}
    return Field(u.space, -u.value)
end

function Base. /(u::Field{S}, a::Number)::Field{S} where {S}
    return Field(u.space, u.value./a)
end

function Base. log(u::Field{S})::Field{S} where {S}
    return Field(u.space, log.(u.value))
end

function Base. exp(u::Field{S})::Field{S} where {S}
    return Field(u.space, exp.(u.value))
end

function reshapeFromTuple(U::NTuple{3, Field})
    return vcat(reshape(U[1]), reshape(U[2]), reshape(U[3]))
end

function reshapeToTuple(space::S, x::Array{T,1})::NTuple{3, Field}  where {S, T}
    U = reshape(x, :, 3)
    return (reshape(space, U[:, 1]), reshape(space, U[:, 2]), reshape(space, U[:, 3]))
end

#--------------------------------------------------------------------
# Auxilliary functions 
#--------------------------------------------------------------------

function F(f::Field{S}, r::Field{S}, ϕ::Field{S})::NTuple{3, Field{S}} where {S}

    # Compute residual at the interior
    resf = DU*DV*log(abs(f)) + (2/r)*(DU*DV*r) + 2*(DU*ϕ)*(DV*ϕ)
    resr = 2*DU*(DV*r) +  (2/r)*(DU*r)*(DV*r) + (f/r)
    resϕ = DU*DV*ϕ + (1/r)*(DU*r)*(DV*ϕ) + (1/r)*(DV*r)*(DU*ϕ)

    # Compute residual at the boundary
    resbndf = B*(f - bndf)  
    resbndr = B*(r - bndr)  
    resbndϕ = B*(ϕ - bndϕ)  

    # Combine the residuals
    finalresf = (I-B)*resf + resbndf
    finalresr = (I-B)*resr + resbndr
    finalresϕ = (I-B)*resϕ + resbndϕ

    return (finalresf, finalresr, finalresϕ)
end

function constraints(f::Field{S}, r::Field{S}, ϕ::Field{S})::NTuple{2, Field{S}} where {S}
    # compute the constraints
    C1 = -(2/r)*DU*DU*r + (2/r)*(1/f)*(DU*r)*(DU*f) - (DU*ϕ)^2
    C2 = -(2/r)*DV*DV*r + (2/r)*(1/f)*(DV*r)*(DV*f) - (DV*ϕ)^2
    return (C1, C2)
end

#--------------------------------------------------------------------
# Schwarzschild Spacetime
#--------------------------------------------------------------------

function find_t_of_UV(U::T, V::T, M::T)::T where {T<:Number}
    @assert V > 0   # ensure you're in region I or II
    @assert U*V < 1 # ensure you don't hit the singularity
    if U*V == 0     # r = 2M 
        t = V       # enforce uniqueness 
    elseif U > 0    # r < 2M
        t = -2M*log(U/V)
    elseif U < 0    # r > 2M
        t = -2M*log(-U/V)
    else
        error("Domain error")
    end
    return t
end

function find_r_of_UV(U::T, V::T, M::T)::T where {T<:Number}
    @assert V > 0       # ensure you're in region I or II
    @assert U*V < 1     # ensure you don't hit the singularity
    if U*V == 0     # r = 2M
        r = 2M
    else            # r < 2M or r > 2M
        f(r) = (r/2M - 1)*exp(r/2M) + U*V
        r    = find_zero(f, 2M)
    end
    @assert r > 0
    return r
end


#--------------------------------------------------------------------
# Generic Spacetime 
#--------------------------------------------------------------------

struct U end
struct V end
S1 = ChebyshevGL{U, 16, Float64}(-8, -6)
S2 = ChebyshevGL{V, 16, Float64}(3, 4)
S  = ProductSpace(S1, S2)
DU, DV = derivative(S)
B = incomingboundary(S)
I = identity(S)

if true
    println("Testing Schwarzschild in vacuum")

    # Construct Initial Data
    M = Float64(1.0)
    r0 = Field(S, (u,v)->v-u)
    ϕ0 = Field(S, (u,v)->exp(-u^2)/0.1 + exp(-v^2)/0.1)
    ru = IDSolver(extractUboundary(r0), extractUboundary(ϕ0))
    rv = IDSolver(extractVboundary(r0), extractVboundary(ϕ0))
    r0 = combineUVboundary(ru, rv)
    f0 = Field(S, (u,v)->1)

    # Boundary conditions
    bndf = B*f0
    bndr = B*r0
    bndϕ = B*ϕ0

    @show bndf

    # Solve
    u = nlsolve(f!, reshapeFromTuple((f0, r0, ϕ0)); autodiff=:forward, show_trace=true, ftol=1e-12)
    fsol, rsol, ϕsol = reshapeToTuple(S, u.zero)
    
    # Check constraints
    C1, C2 = constraints(fsol, rsol, ϕsol)
    @show L2(C1)
    @show L2(C2)
end

