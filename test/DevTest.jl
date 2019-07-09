#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 07-2019
# Solve the wave equation with a fixed background 
# in spherically-symmetric spacetimes 
#--------------------------------------------------------------------

struct U end
struct V end
using LinearAlgebra


function LinearAlgebra. display(A::Operator{S}) where {S}
    display(reshape(A))
end

function LinearAlgebra. display(A::Field{S}) where {S}
    display(A.value)
end

# Construct the space
SU  = ChebyshevGL{U, 2, Float64}(-1, 1)
SV  = ChebyshevGL{V, 2, Float64}(-1, 1)
SUV = ProductSpace(SU, SV)

# Construct the derivatives and the field
DU, DV = derivative(SUV)
f = Field(SUV, (U,V)->1)
r = Field(SUV, (U,V)->U-V)

# Construct the operator
L = DU*DV + ((DV*r)/r)*DU + ((DU*r)/r)*DV

# Constuct the boundary operator
B = incomingboundary(SUV)

# Enforce regularity on axis and boundary conditions at the incoming boundaries
A = axisboundary(SUV) 
M = Operator(A.space, isfinite.(L.value).*L.value)
L = (M + A*(DU*DV)) ⊕ B
# display(L)
println()
@show cond(reshape(L))
@show sort(abs.(eigvals(reshape(L))))

# Now set boundary conditions
b = Field(SUV, (U,V)->sin(U-V)*cos(U+V)) 
ϕ = solve(L, B*b)
display(ϕ)


