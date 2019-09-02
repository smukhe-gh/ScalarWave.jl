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

function extractUboundary(u::Field{ProductSpace{S1,S2}})::Field{S1} where {S1, S2}
    @assert ndims(u.value) == 2
    return Field(u.space.S1, u.value[end, :])
end

function extractVboundary(u::Field{ProductSpace{S1,S2}})::Field{S2} where {S1, S2}
    @assert ndims(u.value) == 2
    return Field(u.space.S2, u.value[:, end])
end

function combineUVboundary(u::Field{S1}, v::Field{S2})::Field{ProductSpace{S1, S2}} where {S1, S2}
    PS = ProductSpace(u.space, v.space)
    w  = Field(PS, (u,v)->rand())    
    w.value[end, :] = u.value
    w.value[:, end] = v.value
    return w
end

#--------------------------------------------------------------------
# Generic Spacetime 
#--------------------------------------------------------------------

using DoubleFloats
struct U end
struct V end

T  = Double64
SU = ChebyshevGL{U, 3, T}(-4, -3)
SV = ChebyshevGL{V, 3, T}( 5,  6)
S  = ProductSpace(SU, SV)

# Schwarzschild Analytic solution
M = T(1.0)
r = Field(S, (u,v)->find_r_of_UV(u,v,M))
f = ((16*M^3)/r)*exp(-r/2M)
ϕ = Field(S, (u,v)->0)

# test constraint equations and it's derivatives with 
# Schwarzschild solution
cu, cv = constraints(f, r, ϕ)
@show L2(cu)
@show L2(cv)

# Now extract the boundaries and check constraints on the line. 
u = Field(S, (u,v)->u)
v = Field(S, (u,v)->v)

@show extractUboundary(f)

exit()

# Now test constraints along lines
clu = lineconstraints(extractUboundary(f),
                      extractUboundary(r),
                      extractUboundary(ϕ))

@show L2(clu)
