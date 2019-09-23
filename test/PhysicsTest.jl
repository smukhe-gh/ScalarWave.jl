#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 10-2019
# Test Physics functions
#--------------------------------------------------------------------

# Schwarzschild spacetime 
# PS = ProductSpace(ChebyshevGL{U, 20, Float64}(-3.0,  -1.0), 
                  # ChebyshevGL{V, 20, Float64}( 1.0,   3.0))

# M  = 1.0
# r0 = Field(PS, (u,v)->find_r_of_UV(u,v,M))
# ϕ0 = Field(PS, (u,v)->0)
# f0 = ((16*M^3)/r0)*exp(-r0/2M)
# a0 = -sqrt(2*f0) 
# n0 = Field(PS, (u,v)->1e-3*rand())

# @show lineconstraint(extractUboundary.((a0, r0, ϕ0), :incoming)...)
# @show lineconstraint(extractVboundary.((a0, r0, ϕ0), :incoming)...)
# @show L2.(constraints(a0, r0, ϕ0))
# @show L2.(constraints(a0 + n0, r0 + n0, ϕ0 + n0))
# (asol, rsol, ϕsol) = compute(computeUboundary((a0, r0, ϕ0)),
                             # computeVboundary((a0, r0, ϕ0)), (a0 + n0, r0 + n0, f0 + n0))
# @show L2(asol - a0)
# @show L2(rsol - r0)
# @show L2(ϕsol - ϕ0)
# println()
# println()
# println()

# Minkowski background + scalar wave
PS = ProductSpace(ChebyshevGL{U, 15, Float64}( 0, 1), 
                  ChebyshevGL{V, 30, Float64}( 2, 3))
r0 = Field(PS, (u,v)->v-u)
ϕ0 = Field(PS, (u,v)->0.6*exp(-(v-4)^2))
f0 = Field(PS, (u,v)->1)
a0 = -sqrt(2*f0) 

# @show lineconstraint(extractUboundary.((a0, r0, ϕ0), :incoming)...)
# @show lineconstraint(extractVboundary.((a0, r0, ϕ0), :incoming)...)
# @show L2.(constraints(a0, r0, ϕ0))
# (asol, rsol, ϕsol) = compute(computeUboundary((a0, r0, ϕ0)),
                             # computeVboundary((a0, r0, ϕ0)), (a0, r0, f0))
# @show L2.(constraints(asol, rsol, ϕsol))

# Minkowski background + scalar wave on Multiple patches
grid = Grid(Float64, (8, 10), (4, 14), (0,1), (2,3))
tree = distribute(grid, (u,v)->0.6*exp(-(v-2.5)^2))

C2A = []
C1A = []
for index in CartesianIndices(tree)
    C1, C2 = L2.(constraints(tree[index]...))
    append!(C1A, C1)
    append!(C2A, C2)
end

@show maximum(C1A)
@show maximum(C2A)

# contourf(extractfield(tree, :ϕ), 10)
# show()
