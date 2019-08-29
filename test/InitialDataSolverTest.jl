#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 07-2019
# Write an initial data solver
# See <https://arxiv.org/pdf/1510.05273.pdf> for the equations
# They do not seem to be consistent with Waugh and Lake 1980 
#--------------------------------------------------------------------

#  Choose f and solve for r
function IDSolverR(f::Field{S}, r::Field{S}, ϕ::Field{S})::Field{S} where {S<:Space{Tag}} where {Tag}

    #  Construct local 1D operators
    D = derivative(ϕ.space)
    I = identity(ϕ.space)

    # Specify boundary values at the end points for r
    B = incomingboundary(ϕ.space) + outgoingboundary(ϕ.space)
    b = B*r 

    # Construct operator for r after choosing the gauge for f
    # and solve the 1D equation
    A = 2*(D*D - (1/f)*(D*f)*(D)) + ((D*ϕ)^2)*I
    rsolved = solve(A ⊕ B, b)

    return rsolved
end

function IDSolverF(f::Field{S}, r::Field{S}, ϕ::Field{S}, symbol::Symbol)::Field{S} where {S<:Space{Tag}} where {Tag}

    #  Construct local 1D operators
    D = derivative(ϕ.space)
    I = identity(ϕ.space)

    # Construct RHS to enforce boundary conditions
    if symbol==:U   # check the consistency 
        B = incomingboundary(ϕ.space)
    else
        B = outgoingboundary(ϕ.space)
    end
    b  = B*f 

    # Construct operator for f and solve the 1D equation
    A = 2*((D*D*r)*I - (D*r)*D) + (r*(D*ϕ)^2)*I
    fsolved = solve(A ⊕ B, b)

    return fsolved
end

function extractUboundary(u::Field{ProductSpace{S1,S2}})::Field{S1} where {S1, S2}
    return Field(u.space.S1, u.value[1,:])
end

function extractVboundary(u::Field{ProductSpace{S1,S2}})::Field{S2} where {S1, S2}
    return Field(u.space.S2, u.value[:,1])
end

function combineUVboundary(u::Field{S1}, v::Field{S2})::Field{ProductSpace{S1, S2}} where {S1, S2}
    PS = ProductSpace(u.space, v.space)
    w  = Field(PS, (u,v)->1)    # This is a strange choice, but should not affect our solve as we don't enforce it on the interior. 
    w.value[1,:] = u.value
    w.value[:,1] = v.value
    return w
end

function F(f::Field{S}, r::Field{S}, ϕ::Field{S})::NTuple{3, Field{S}} where {S}
    δf = DU*DV*log(abs(f)) + (2/r)*(DU*DV*r) + 2*(DU*ϕ)*(DV*ϕ)
    δr = 2*DU*(DV*r) +  (2/r)*(DU*r)*(DV*r) + (f/r)
    δϕ = DU*DV*ϕ + (1/r)*(DU*r)*(DV*ϕ) + (1/r)*(DV*r)*(DU*ϕ)
    return (δf, δr, δϕ)
end

function C(f::Field{S}, r::Field{S}, ϕ::Field{S})::NTuple{2, Field{S}} where {S}
    δu = 2*(DU*DU*r - (1/f)*(DU*f)*(DU*r)) + r*(DU*ϕ)^2
    δv = 2*(DV*DV*r - (1/f)*(DV*f)*(DV*r)) + r*(DV*ϕ)^2
    return (δu, δv)
end

function lineC(f::Field{S}, r::Field{S}, ϕ::Field{S})::Field{S} where {S}
    D = derivative(ϕ.space)
    return 2*(D*D*r - (1/f)*(D*f)*(D*r)) + r*(D*ϕ)^2
end

#--------------------------------------------------------------------
# Generic Spacetime 
#--------------------------------------------------------------------

struct U end
struct V end

T  = Complex
S1 = ChebyshevGL{U, 16, T}(3, 4)
S2 = ChebyshevGL{V, 16, T}(3, 4)

S  = ProductSpace(S1, S2)
DU, DV = derivative(S)
B = incomingboundary(S)
I = identity(S)

# Check constraint on incoming conditions
ru0 = Field(S1, u->exp(im*u/sqrt(2)))
fu0 = Field(S1, u->1)
ϕu0 = Field(S1, u->u)
@show L2(lineC(fu0, ru0, ϕu0))

# Now promote them to 2D
ru02D = Field(S, (u,v)->exp(im*(u+v)/sqrt(2)))
fu02D = Field(S, (u,v)->1)
ϕu02D = Field(S, (u,v)->u+v)
@show L2.(C(fu02D, ru02D, ϕu02D))

# Recover the analytic solutions using the solver
rsolved = IDSolverR(fu0, ru0, ϕu0) 
@test rsolved ≈ ru0

