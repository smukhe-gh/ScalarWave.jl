#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 10-2019
#--------------------------------------------------------------------

PS = ProductSpace(ChebyshevGL{U, 18, Float64}(0, 1), 
                  ChebyshevGL{V, 18, Float64}(0, 1))

p0 = 0.01
a0 = Field(PS, (u,v)->1) 
η0 = Field(PS, (u,v)->(v-u)/2) 
ϕ0 = Field(PS, (u,v)->p0*exp(-(v-1/2)^2) + p0*exp(-(u-1/2)^2))

(abnd, ηbnd, ϕbnd) = combineUVboundary.(computeUboundary((a0, η0, ϕ0)), 
                                        computeVboundary((a0, η0, ϕ0)), :incoming)

@show lineconstraint(extractUboundary.((abnd, ηbnd, ϕbnd), :incoming)...)
@show lineconstraint(extractVboundary.((abnd, ηbnd, ϕbnd), :incoming)...)


