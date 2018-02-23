#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

struct Patch
    loc::Array{Int,1}
    value::Array{Float64,2}
end

struct Boundary
    kind::Int
    value::Array{Float64,1}
end

struct LocalC
    loc::Array{Float64,1}
end

struct GlobalC
    loc::Array{Float64,1}
end

