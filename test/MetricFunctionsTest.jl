#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 04-2018
# Test metric functions 
#--------------------------------------------------------------------

struct U end
struct V end
struct UV end

#--------------------------------------------------------------------
# Define boundary and the product space
#--------------------------------------------------------------------
P1, P2 = 2, 2
M   = 1.0
SUV = ProductSpace{GaussLobatto{U,P1}, GaussLobatto{V,P2}}

#--------------------------------------------------------------------
# Define metric functions 
#--------------------------------------------------------------------

𝒈tt = Field(SUV, (u,v)->1)  
𝒈rr = Field(SUV, (u,v)->2) 
𝒈θθ = Field(SUV, (u,v)->3) 
𝒈ϕϕ = Field(SUV, (u,v)->4) 
𝒈rθ = Field(SUV, (u,v)->5) 
𝒈rϕ = Field(SUV, (u,v)->6) 
𝒈tr = Field(SUV, (u,v)->7)  
𝒈tθ = Field(SUV, (u,v)->8) 
𝒈tϕ = Field(SUV, (u,v)->9) 
𝒈θϕ = Field(SUV, (u,v)->10) 

𝕘 = Metric{dd, 4}([𝒈tt, 𝒈tr, 𝒈tθ, 𝒈tϕ,
                        𝒈rr, 𝒈rθ, 𝒈rϕ, 
                             𝒈θθ, 𝒈θϕ, 
                                  𝒈ϕϕ])

@show typeof(𝕘)
@show dim(𝕘)
@test 𝕘[2,3] ==  𝕘[3,2]
@show metricdet(𝕘)

