#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2018
# Define concrete datatypes for GR
#--------------------------------------------------------------------

struct Derivative{Tag, D}
    components::Array{T, 1} where {T}
end

struct Metric{Tag, D}
    components::Array{T, 1} where {T}
end

mutable struct Christoffel{Tag, D}
    components::Array{Field, 3}
end

mutable struct Ricci{Tag, D}
    components::Array{Field, 1}
end
