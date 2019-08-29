#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 07-2019
# Write a non-linear solver from scratch for Einstein's field
# equations
#--------------------------------------------------------------------

using DualNumbers
using NLsolve
using Roots

#--------------------------------------------------------------------
# Core functions 
#--------------------------------------------------------------------

struct Residual{S}
    f::Field{S}
    r::Field{S}
    ϕ::Field{S}
end

struct State{S}
    f::Field{S}
    r::Field{S}
    ϕ::Field{S}
end

function f!(fvec::Array{T,1}, x::Array{T,1}) where {T}
    fvec[:] = reshape(F(reshape(x, S)))
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

function Base. reshape(x::Array{T,1}, space::S)::State{S} where {T, S <: ProductSpace}
    xshape = reshape(x, (:, 3))
    return State(Field(space, xshape[:, 1]), 
                 Field(space, xshape[:, 2]),
                 Field(space, xshape[:, 3]))
end

function Base. reshape(R::Residual{S}) where {S}
    return vcat(reshape(R.f), reshape(R.ϕ), reshape(R.ϕ))
end

function Base. reshape(U::State{S}) where {S}
    return vcat(reshape(U.f), reshape(U.ϕ), reshape(U.ϕ))
end

function Base. +(R::Residual{S}, Q::Residual{S})::Residual{S} where {S}
    return Residual(R.f + Q.f, R.r + Q.r, R.ϕ + Q.ϕ)
end

function Base. *(A::Operator{S}, R::Residual{S})::Residual{S} where {S}
    return Residual(A*R.f, A*R.r, A*R.ϕ) 
end

#--------------------------------------------------------------------
# Auxilliary functions 
#--------------------------------------------------------------------

function F(state::State{S})::Residual{S} where {S}
    (f, r, ϕ)  = (state.f, state.r, state.ϕ)

    resf = DU*DV*log(abs(f)) + (2/r)*(DU*DV*r) + 2*(DU*ϕ)*(DV*ϕ)
    resr = 2*DU*(DV*r) +  (2/r)*(DU*r)*(DV*r) + (f/r)
    resϕ = DU*DV*ϕ + (1/r)*(DU*r)*(DV*ϕ) + (1/r)*(DV*r)*(DU*ϕ)

    InteriorResidual = Residual(resf, resr, resϕ)
    BoundaryResidual = Residual(f - bndf, r - bndr, ϕ  - bndϕ)

    return (I-B)*InteriorResidual + B*BoundaryResidual
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
# Test solve
#--------------------------------------------------------------------

struct Q end
S1 = ChebyshevGL{Q, 3, Float64}(-8, -4)
S2 = ChebyshevGL{Q, 3, Float64}(3, 4)
S  = ProductSpace(S1, S2)
DU, DV = derivative(S)
B = incomingboundary(S)
I = identity(S)

if true
    println("Testing Schwarzschild in vacuum")

    # Analytic solution
    M = Float64(1.0)
    r = Field(S, (u,v)->find_r_of_UV(u,v,M))
    f = ((16*M^3)/r)*exp(-r/2M)
    ϕ = Field(S, (u,v)->0)

    # Boundary conditions
    bndf = f
    bndr = r
    bndϕ = ϕ
    
    residual = F(State(f, r, ϕ))
    @show L2(residual.f)
    @show L2(residual.r)
    @show L2(residual.ϕ)

    # Iniital guess
    fac = 1e-1
    noise = Field(S, (u,v)->fac*rand())
    @show maximum(noise)
    
    u = nlsolve(f!, reshape(State(f + noise, r + noise, ϕ + noise)); autodiff=:forward, show_trace=true, ftol=1e-9)
    solvedState = reshape(u.zero, S)
    @show L2(solvedState.f-f)
    @show L2(solvedState.r-r)
    @show L2(solvedState.ϕ-ϕ)
end
