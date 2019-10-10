#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 07-2019
# Construct different residual functions to check what works
# on axis
#--------------------------------------------------------------------

# Scaled grr = η^2 ω and store ω on grid
function F(a::Field{S}, ω::Field{S}, ϕ::Field{S})::NTuple{3, Field{S}} where {S}
    F1 = (η*ω)*(DU*(DV*ϕ)) + (DU*ϕ)*(η*(DV*ω) + ω*(DV*η)) + (DV*ϕ)*(η*(DU*ω) + ω*(DU*η))
    F2 = (ω^2)*(η*(DU*(DV*η)) + (DU*η)*(DV*η)) + (η^2)*(DU*ω)*(DV*ω) + (ω*η)*(η*(DU*(DV*ω)) + 2*(DU*η)*(DV*ω) + 2*(DU*ω)*(DV*η)) + (1/4)*(a^2)
    F3 = (1/a)*(DU*(DV*a)) - (1/a^2)*(DU*a)*(DV*a) + (1/ω)*(DU*(DV*ω)) + 4pi*(DU*ϕ)*(DV*ϕ)
    return F1, F2, F3
end

# Residual functions with regularity condition on the axis 
function F(a::Field{S}, r::Field{S}, ϕ::Field{S})::NTuple{3, Field{S}} where {S}
    F1 = r*(DU*(DV*ϕ)) + (DU*r)*(DV*ϕ) + (DV*r)*(DU*ϕ)
    F2 = r*(DU*(DV*r)) + (DU*r)*(DV*r) + (1/4)*(a^2)
    F3 = (1/a)*(DU*(DV*a)) - (1/a^2)*(DU*a)*(DV*a) + (1/r)*(DU*(DV*r)) + 4pi*(DU*ϕ)*(DV*ϕ)
    F3onAxis = (1/a)*(DU*(DV*a)) - (1/a^2)*(DU*a)*(DV*a) + 4pi*(DU*ϕ)*(DV*ϕ)

    for index in CartesianIndices(F3.value) 
        if index.I[1] == index.I[2]
            F3.value[index] = F3onAxis.value[index]
        end
    endh

    return (F1, F2, F3)
end

# Standard residual functions
function F(a::Field{S}, r::Field{S}, ϕ::Field{S})::NTuple{3, Field{S}} where {S}
    F1 = r*(DU*(DV*ϕ)) + (DU*r)*(DV*ϕ) + (DV*r)*(DU*ϕ)
    F2 = r*(DU*(DV*r)) + (DU*r)*(DV*r) + (1/4)*(a^2)
    F3 = (1/a)*(DU*(DV*a)) - (1/a^2)*(DU*a)*(DV*a) + (1/r)*(DU*(DV*r)) + 4pi*(DU*ϕ)*(DV*ϕ)
    return (F1, F2, F3)
end

# Residual for Minkowski with singular term deleted
function F(a::Field{S}, r::Field{S}, ϕ::Field{S})::NTuple{3, Field{S}} where {S}
    DU, DV = derivative(a.space)
    F1 = 0*r
    F2 = r*(DU*(DV*r)) + (DU*r)*(DV*r) + (1/4)*(a^2)
    F3 = (1/a)*(DU*(DV*a)) - (1/a^2)*(DU*a)*(DV*a)
    return (F1, F2, F3)
end
