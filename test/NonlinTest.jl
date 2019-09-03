#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 07-2019
# Write a non-linear solver from scratch for Einstein's field
# equations
#--------------------------------------------------------------------

using DualNumbers
using NLsolve
using Roots
using DoubleFloats

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

#--------------------------------------------------------------------
# Auxilliary functions 
#--------------------------------------------------------------------

function Base. log(u::Field{S})::Field{S} where {S}
    return Field(u.space, log.(u.value))
end

function Base. exp(u::Field{S})::Field{S} where {S}
    return Field(u.space, exp.(u.value))
end

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

function C(f::Field{S}, r::Field{S}, ϕ::Field{S})::NTuple{2, Field{S}} where {S}
    CU = 2*((1/f)*(DU*f)*(DU*r) - DU*DU*r)/r
    CV = 2*((1/f)*(DV*f)*(DV*r) - DV*DV*r)/r
    return (CU, CV)
end

function reshapeFromTuple(U::NTuple{3, Field})
    return vcat(reshape(U[1]), reshape(U[2]), reshape(U[3]))
end

function reshapeToTuple(space::S, x::Array{T,1})::NTuple{3, Field}  where {S, T}
    U = reshape(x, :, 3)
    return (reshape(space, U[:, 1]), reshape(space, U[:, 2]), reshape(space, U[:, 3]))
end

#--------------------------------------------------------------------
# Test solve
#--------------------------------------------------------------------

struct Q end
T  = Double64
S1 = ChebyshevGL{Q, 18, T}(-8, -4)
S2 = ChebyshevGL{Q, 13, T}( 3,  5)
S  = ProductSpace(S1, S2)
DU, DV = derivative(S)
B = outgoingboundary(S)
I = identity(S)

#--------------------------------------------------------------------
# Schwarzschild Spacetime
#--------------------------------------------------------------------

if true
    println("Testing Schwarzschild in vacuum")

    # Analytic solution
    M = T(1.0)
    r = Field(S, (u,v)->find_r_of_UV(u,v,M))
    f = ((16*M^3)/r)*exp(-r/2M)
    ϕ = Field(S, (u,v)->0)

    # Boundary conditions
    bndf = B*f
    bndr = B*r
    bndϕ = B*ϕ
    
    resf, resr, resϕ = F(f, r, ϕ)
    CU, CV = C(f, r, ϕ)
    @show L2(CU)
    @show L2(CV)

    # Iniital guess
    fac = 1e-1
    noise = Field(S, (u,v)->fac*rand())
    @show maximum(noise)
    
    u = nlsolve(f!, reshapeFromTuple((f + noise, r + noise, ϕ + noise)); autodiff=:forward, show_trace=true, ftol=1e-12)
    fsol, rsol, ϕsol = reshapeToTuple(S, u.zero)
    @show L1(fsol-f)
    @show L1(rsol-r)
    @show L1(ϕsol-ϕ)

    CU, CV = C(fsol, rsol, ϕsol)
    @show L2(CU)
    @show L2(CV)
end


