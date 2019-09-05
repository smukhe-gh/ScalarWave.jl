#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 07-2019
# Write a non-linear solver from scratch for Einstein's field
# equations
#--------------------------------------------------------------------

using NLsolve
using PyPlot

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
    DU, DV = derivative(ϕ.space)
    δu = 2*(DU*DU*r - (1/f)*(DU*f)*(DU*r)) + r*(DU*ϕ)^2
    δv = 2*(DV*DV*r - (1/f)*(DV*f)*(DV*r)) + r*(DV*ϕ)^2
    return (δu, δv)
end

function lineconstraints(f::Field{S}, r::Field{S}, ϕ::Field{S})::Field{S} where {S}
    D = derivative(ϕ.space)
    return 2*(D*D*r - (1/f)*(D*f)*(D*r)) + r*(D*ϕ)^2
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
S1 = ChebyshevGL{U, 23, T}(-8, -6)
S2 = ChebyshevGL{V, 23, T}(3,   5)
S  = ProductSpace(S1, S2)
DU, DV = derivative(S)
B = outgoingboundary(S)
I = identity(S)

# Generic spacetime
M = T(1.0)
f = Field(S, (u,v)->1)
r = Field(S, (u,v)->v-u)
ϕ = Field(S, (u,v)->exp(-(v-4)^2))  # Incoming wave travelling in the u-direction

# Schwarzschild spacetime 
# M = T(1.0)
# r = Field(S, (u,v)->find_r_of_UV(u,v,M))
# f = ((16*M^3)/r)*exp(-r/2M)
# ϕ = Field(S, (u,v)->0)

# Construct constraint-satisfying initial data by fixing f and solving for r
rsolvedU = IDSolverR(extractUboundary(f), extractUboundary(r), extractUboundary(ϕ))
rsolvedV = IDSolverR(extractVboundary(f), extractVboundary(r), extractVboundary(ϕ))
rsolved  = combineUVboundary(rsolvedU, rsolvedV)

@show L2(lineconstraints(extractUboundary(f), extractUboundary(rsolved), extractUboundary(ϕ)))
@show L2(lineconstraints(extractVboundary(f), extractVboundary(rsolved), extractVboundary(ϕ)))

# Construct constraint-satisfying initial data by fixing r and solving for f
# fsolvedU = IDSolverF(extractUboundary(f), extractUboundary(r), extractUboundary(ϕ))
# fsolvedV = IDSolverF(extractVboundary(f), extractVboundary(r), extractVboundary(ϕ))
# fsolved  = combineUVboundary(fsolvedU, fsolvedV)

# @show L2(lineconstraints(extractUboundary(fsolved), extractUboundary(r), extractUboundary(ϕ)))
# @show L2(lineconstraints(extractVboundary(fsolved), extractVboundary(r), extractVboundary(ϕ)))


# Derive boundary conditions from initial data
bndf = B*fguess
bndr = B*rsolved
bndϕ = B*ϕguess

# Compute constraints before solve
CU, CV = constraints(fguees, rsolved, ϕguess)
@show L2(CU)
@show L2(CV)

#--------------------------------------------------------------------
# Now start the nonlinear solve
#--------------------------------------------------------------------
u = nlsolve(f!, reshapeFromTuple((fguess, rguess, ϕguess)); autodiff=:forward, show_trace=true, ftol=1e-9)
fsol, rsol, ϕsol = reshapeToTuple(S, u.zero)

CU, CV = constraints(fsol, rsol, ϕsol)
@show L2(CU)
@show L2(CV)

# Check constraint violations at the incoming boundary
# Note that residuals here are of the order of 1e-11
@show L2(lineconstraints(extractUboundary(fsol), extractUboundary(rsol), extractUboundary(ϕsol)))
@show L2(lineconstraints(extractVboundary(fsol), extractVboundary(rsol), extractVboundary(ϕsol)))

#--------------------------------------------------------------------
# Tests
#--------------------------------------------------------------------
# [1] Plot r coordinate
#       -- The r coordinate is not monotonic. This might not be physical. 
#       -- The solver seems to be highly sensitive towards the initial guess. Thus we might need to combine U and V boundaries with more care.
#          > Solve for f specifying r. 
# [2] Test residual at the incoming boundaries
#       -- FIXME: Residual is 1e-11 at the boundaries. Constraint violations however are large at the boundaries.
#          > Check constraint violations at the boundaries again. How can residual be small and the constraints violations be large? 
# [3] Plot the solutions after solve [Done]
# [4] Plot residual after solve [Done--not very informative.]
# [5] Plot constraint violation after solve [Done]
# [6] Understand where the constraint violations are the largest--experiment with size of patch and number of points. What to make of it? 
#       -- The scalar wave should contract under self gravity. In my case, it seems to be diverging. Need more physical intuition. 


