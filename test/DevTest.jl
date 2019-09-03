#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 07-2019
# Write a non-linear solver from scratch for Einstein's field
# equations
#--------------------------------------------------------------------

using DualNumbers
using NLsolve
using Roots
using PyPlot
import DoubleFloats.Double64

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

function reshapeFromTuple(U::NTuple{3, Field})
    return vcat(reshape(U[1]), reshape(U[2]), reshape(U[3]))
end

function reshapeToTuple(space::S, x::Array{T,1})::NTuple{3, Field}  where {S, T}
    U = reshape(x, :, 3)
    return (reshape(space, U[:, 1]), reshape(space, U[:, 2]), reshape(space, U[:, 3]))
end

#--------------------------------------------------------------------
# Initial data solver
#--------------------------------------------------------------------

function constraints(f::Field{S}, r::Field{S}, ϕ::Field{S})::NTuple{2, Field{S}} where {S}
    DU, DV = derivative(ϕ.space)
    δu = 2*(DU*DU*r - (1/f)*(DU*f)*(DU*r)) + r*(DU*ϕ)^2
    δv = 2*(DV*DV*r - (1/f)*(DV*f)*(DV*r)) + r*(DV*ϕ)^2
    return (δu, δv)
end

function lineconstraints(f::Field{S}, r::Field{S}, ϕ::Field{S})::Field{S} where {S}
    D = derivative(ϕ.space)
    return 2*(D*D*r - (1/f)*(D*f)*(D*r)) + r*(D*ϕ)^2
end

function extractUboundary(u::Field{ProductSpace{S1,S2}})::Field{S2} where {S1, S2}
    @assert ndims(u.value) == 2
    return Field(u.space.S2, u.value[end, :])
end

function extractVboundary(u::Field{ProductSpace{S1,S2}})::Field{S1} where {S1, S2}
    @assert ndims(u.value) == 2
    return Field(u.space.S1, u.value[:, end])
end

function combineUVboundary(ubnd::Field{SV}, vbnd::Field{SU})::Field{ProductSpace{SU, SV}} where {SU, SV}
    PS = ProductSpace(vbnd.space, ubnd.space)
    w  = Field(PS, (u,v)->rand())    
    w.value[end, :] = ubnd.value
    w.value[:, end] = vbnd.value
    return w
end

#  Choose f and solve for r
function IDSolverR(f::Field{S}, r::Field{S}, ϕ::Field{S})::Field{S} where {S<:Space{Tag}} where {Tag}
    #  Construct local 1D operators
    D = derivative(ϕ.space)
    I = identity(ϕ.space)
    A = 2*(D*D - (1/f)*(D*f)*(D)) + ((D*ϕ)^2)*I

    # Specify boundary values at the end points for r
    B = incomingboundary(ϕ.space) + outgoingboundary(ϕ.space)
    b = B*r 

    # Construct operator for r after choosing the gauge for f
    # and solve the 1D equation
    rsolved = solve(A ⊕ B, b)
    return rsolved
end

# Choose r and solve for f
function IDSolverF(f::Field{S}, r::Field{S}, ϕ::Field{S})::Field{S} where {S<:Space{Tag}} where {Tag}
    #  Construct local 1D operators
    D = derivative(ϕ.space)
    I = identity(ϕ.space)
    A = 2*((D*D*r)*I - (D*r)*D) + (r*(D*ϕ)^2)*I

    # Construct RHS to enforce boundary conditions
    B = outgoingboundary(ϕ.space)
    b = B*f 

    # Construct operator for f and solve the 1D equation
    fsolved = solve(A ⊕ B, b)
    return fsolved
end

#--------------------------------------------------------------------
# Solve for an arbitrary spacetime
#--------------------------------------------------------------------

struct U end
struct V end
T  = Float64 #Double64
S1 = ChebyshevGL{U, 33, T}(-8, -4)
S2 = ChebyshevGL{V, 33, T}( 3,  5)
S  = ProductSpace(S1, S2)
DU, DV = derivative(S)
B = outgoingboundary(S)
I = identity(S)

# Generic spacetime
M = T(1.0)
r = Field(S, (u,v)->v-u)
f = Field(S, (u,v)->1)
ϕ = Field(S, (u,v)->exp(-(v-4)^2))  # Incoming wave travelling in the u-direction

# Schwarzschild spacetime 
# M = T(1.0)
# r = Field(S, (u,v)->find_r_of_UV(u,v,M))
# f = ((16*M^3)/r)*exp(-r/2M)
# ϕ = Field(S, (u,v)->0)

# Construct constraint-satisfying initial data
rsolvedU = IDSolverR(extractUboundary(f), 
                     extractUboundary(r), 
                     extractUboundary(ϕ))

# plot(extractUboundary(r))
# plot(rsolvedU)
# show()

rsolvedV = IDSolverR(extractVboundary(f), 
                     extractVboundary(r), 
                     extractVboundary(ϕ))

# plot(extractVboundary(r))
# plot(rsolvedV)
# show()

rsolved = combineUVboundary(rsolvedU, rsolvedV)

@show L2(lineconstraints(extractUboundary(f), rsolvedU, extractUboundary(ϕ)))
@show L2(lineconstraints(extractVboundary(f), rsolvedV, extractVboundary(ϕ)))

# Derive boundary conditions from initial data
bndf = B*f
bndr = B*rsolved
bndϕ = B*ϕ

# Compute constraints before solve
resf, resr, resϕ = F(f, rsolved, ϕ)
CU, CV = constraints(f, rsolved, ϕ)
@show L2(CU)
@show L2(CV)

# Now start the nonlinear solve
u = nlsolve(f!, reshapeFromTuple((f, r, ϕ)); autodiff=:forward, show_trace=true, ftol=1e-10)
fsol, rsol, ϕsol = reshapeToTuple(S, u.zero)

CU, CV = constraints(fsol, rsol, ϕsol)
@show L2(CU)
@show L2(CV)



#--------------------------------------------------------------------
# Tests
#--------------------------------------------------------------------
# [1] Plot r coordinate
#       -- The r coordinate is not monotonic. This might not be physical. 
#       -- The solver seems to be highly sensitive towards the initial guess. Thus we might need to combine U and V boundaries with more care.
#          > Solve for f specifying r. 
# [2] Test residual at the incoming boundaries
#       -- Residual is 1e-11 at the boundaries. Constraint violations however are large at the boundaries.
#          > Check constraint violations at the boundaries again. How can residual be small and the constraints violations be large? 
# [3] Plot residual after solve
# [4] Plot constraint violation after solve
# [5] Understand where the constraint violations are the largest--experiment with size of patch and number of points 


