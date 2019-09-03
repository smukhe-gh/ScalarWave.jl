#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 07-2019
# Write an initial data solver
# See <https://arxiv.org/pdf/1510.05273.pdf> for the equations
# They do not seem to be consistent with Waugh and Lake 1980 
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

#--------------------------------------------------------------------
# Generic Spacetime 
#--------------------------------------------------------------------

import DoubleFloats.Double64
struct U end
struct V end

T  = Double64
SU = ChebyshevGL{U, 6, T}(-4, -3) 
SV = ChebyshevGL{V, 6, T}( 5,  6)
S  = ProductSpace(SU, SV)

# Schwarzschild Analytic solution
M = T(1.0)
r = Field(S, (u,v)->find_r_of_UV(u,v,M))
f = ((16*M^3)/r)*exp(-r/2M)
ϕ = Field(S, (u,v)->0)

# test constraint equations and it's derivatives with 
# Schwarzschild solution
# cu, cv = constraints(f, r, ϕ)
# @show L2(cu)
# @show L2(cv)

# Now test constraints along lines
# clu = lineconstraints(extractUboundary(f),
                      # extractUboundary(r),
                      # extractUboundary(ϕ))

# clv = lineconstraints(extractVboundary(f),
                      # extractVboundary(r),
                      # extractVboundary(ϕ))
# @show L2(clu)
# @show L2(clv)

#--------------------------------------------------------------------
# Now test the solvers
#--------------------------------------------------------------------

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

rsolvedU = IDSolverR(extractUboundary(f), 
                     extractUboundary(r), 
                     extractUboundary(ϕ))

rsolvedV = IDSolverR(extractVboundary(f), 
                     extractVboundary(r), 
                     extractVboundary(ϕ))

# rsolved = combineUVboundary(rsolvedU, rsolvedV)
# display(outgoingboundary(S)*(rsolved - r))

# @show L2(extractUboundary(r) - rsolved)

# fsolved = IDSolverF(extractVboundary(f), 
                    # extractVboundary(r), 
                    # extractVboundary(ϕ))
# @show L2(extractVboundary(f) - fsolved)

#--------------------------------------------------------------------
# Now test with generic spacetimes
#--------------------------------------------------------------------

M = T(1.0)
r = Field(S, (u,v)->u-v)
f = Field(S, (u,v)->1)
ϕ = Field(S, (u,v)->exp(-v^2) + exp(-u^2))

rsolvedU = IDSolverR(extractUboundary(f), 
                     extractUboundary(r), 
                     extractUboundary(ϕ))

rsolvedV = IDSolverR(extractVboundary(f), 
                     extractVboundary(r), 
                     extractVboundary(ϕ))

@show L2(lineconstraints(extractUboundary(f), rsolvedU, extractUboundary(ϕ)))
@show L2(lineconstraints(extractVboundary(f), rsolvedV, extractVboundary(ϕ)))

@show rsolvedU.value
@show rsolvedV.value
rsolved = combineUVboundary(rsolvedU, rsolvedV)
display(rsolved)


