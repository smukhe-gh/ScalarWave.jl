#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 07-2019
# Simulate the collapse of a field
#--------------------------------------------------------------------

struct U end
struct V end

PS = ProductSpace(ChebyshevGL{U, 15, Float64}(0, 1), 
                  ChebyshevGL{V, 60, Float64}(2, 3))

a0, r0, ϕ0 = initialguess(PS, (u,v)->0.6*(exp(-v-2.5)^2))
bnda, bndr, bndϕ = initialdatasolver(a0, r0, ϕ0)
a, r, ϕ = nonlinearsolver(PS, (bnda, bndr, bndϕ), (a0, r0, ϕ0))

