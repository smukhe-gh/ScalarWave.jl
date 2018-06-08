#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 04-2018
# Types for Schwarzschild geometry
#--------------------------------------------------------------------

@enum Direction du=1 dv=2

struct Params
    p::Tuple
    size::Tuple
    xmin::Tuple
    xmax::Tuple
    mass::Float64
end

struct Grid
    params::Params
    U::Array{Float64,2}
    V::Array{Float64,2}
    t::Array{Float64,2}
    r::Array{Float64,2}
    drdU::Array{Float64,2}
    drdV::Array{Float64,2}
    dtdU::Array{Float64,2}
    dtdV::Array{Float64,2}
    ddrdUdU::Array{Float64,2}
    ddrdUdV::Array{Float64,2}
    ddrdVdV::Array{Float64,2}
    ddtdUdU::Array{Float64,2}
    ddtdUdV::Array{Float64,2}
    ddtdVdV::Array{Float64,2}
end

struct VarList 
    g01::Array{Float64,2}
    g10::Array{Float64,2}
    g22::Array{Float64,2}
    g33::Array{Float64,2}
    detg::Array{Float64,2}
    scalar::Array{Float64,2}
end
