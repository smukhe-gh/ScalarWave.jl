#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2018
# Distorted Minkowski on multiple patches
#--------------------------------------------------------------------

struct CoordinateRange{T<:Real}
    min::T
    max::T
end

# XXX: Should you keep the patches to be mutable? 
mutable struct Patch
    U::CoordinateRange
    V::CoordinateRange
    value::Field
end

function distribute{T<:Integer}(fbndr::Function, fbndc::Function, frhs::Function, Nx::T, Ny::T, M::T)::Dict{Array{Int,1}, Patch}
    dbase  = Dict{Array{Int, 1}, Patch}()
    for i in 2:2M, k in i-min(i-1,M):min(i-1,M)
        loc  = [k, i-k]
        rhs  = RHS(frhs, Nx, Ny, M, loc)
        bndx = (loc[2]==1) ? (getPatchIC(fbndr, 0, Nx, M, loc[1])) : (getPatchBnd(dbase[loc-[0,1]], 0))
        bndy = (loc[1]==1) ? (getPatchIC(fbndc, 1, Ny, M, loc[2])) : (getPatchBnd(dbase[loc-[1,0]], 1))
        dbase[loc] = calcPatch(bndx, bndy, rhs, dop, bop, loc)
    end
    return dbase
end

#--------------------------------------------------------------------
# Define global quanties [or quantities local to all patches]
#--------------------------------------------------------------------

P1, P2 = 20, 20
SUV = ProductSpace{GaussLobatto{U,P1}, GaussLobatto{V,P2}}
𝔹 = boundary(Null, SUV)

#--------------------------------------------------------------------
# Define patch local quanties
#--------------------------------------------------------------------

# Define coordinate range local to a patch
U = CoordinateRange(1, -1)
V = CoordinateRange(1, -1)

# Now compute coordinates local to a patch
u = Field(SUV, (u,v)->((U.max - U.min)/2)*u + ((U.max + U.min)/2))
v = Field(SUV, (u,v)->((V.max - V.min)/2)*v + ((V.max + V.min)/2))
Ω = Field(SUV, (u,v)->(pi/8)*cospi(u/2)*cospi(v/2))
𝒖 =  u*cos(Ω) + v*sin(Ω)
𝒗 = -u*sin(Ω) + v*cos(Ω)

𝔻𝒗, 𝔻𝒖 = derivativetransform(SUV, 𝒖, 𝒗)

#--------------------------------------------------------------------
# Set boundary conditions [if boundary patch]
#--------------------------------------------------------------------
ρ = 0 
𝕤 = exp(-((𝒖^2)/0.1)) 
𝕓 = 𝔹*𝕤

#--------------------------------------------------------------------
# Construct the wave operator in curved spacetime
#--------------------------------------------------------------------
guu = Field(SUV, (u,v)-> 0)
guv = Field(SUV, (u,v)->-2)
gvv = Field(SUV, (u,v)-> 0)

(𝕘𝒖𝒖, 𝕘𝒖𝒗, 𝕘𝒗𝒗) = inversemetrictransform(guu, guv, gvv, 𝒖, 𝒗) 
𝕘   = [𝕘𝒖𝒖 𝕘𝒖𝒗; 𝕘𝒖𝒗 𝕘𝒗𝒗]
𝔻   = [𝔻𝒖, 𝔻𝒗]
𝕃   = 𝕘𝒖𝒗*𝔻𝒖*𝔻𝒗 + 𝕘𝒖𝒗*𝔻𝒗*𝔻𝒖

#--------------------------------------------------------------------
# Solve the system [also check the condition number and eigen values]
#--------------------------------------------------------------------
𝕨 = solve(𝕃 + 𝔹, ρ + 𝕓) 
drawpatch(𝕨, "../output/minkowski-distorted")
