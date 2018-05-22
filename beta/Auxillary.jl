#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 04-2018
# Compactified coordinates Metric functions
#--------------------------------------------------------------------

struct Params
    size::Tuple
    ulim::Tuple
    vlim::Tuple
    mass::Float64
    theta::Float64
end

struct Grid{T<:Union{Array{Float64,2}, Float64}}
    params::Params
    t_of_UV::T
    r_of_UV::T
    drdU::T
    drdV::T
    ddrdUdV::T
    ddrdUdU::T
    ddrdVdV::T
end

struct Metric{T<:Union{Array{Float64,2}, Float64}}
    g01::T
    g10::T
    g22::T
    g33::T
    detg::T
    riccis::T
    R11::T
    R22::T
    R33::T
    R41::T
    R44::T
end

function convertSoA2AoS(params::Params, SoAmetric::Metric)::Array{Metric,2}
    AoSmetric = Array{Metric}(params.size) 
    for index in CartesianRange(params.size)
        AoSmetric[index] = Metric(SoAmetric.g01[index], 
                                  SoAmetric.g10[index], 
                                  SoAmetric.g22[index], 
                                  SoAmetric.g33[index], 
                                  SoAmetric.detg[index], 
                                  SoAmetric.riccis[index]) 
    end
    return AoSmetric
end

function convertAoS2SoA(params::Params, AoSmetric::Array{ScalarWave.Metric,4})::Metric 
    SoAmetric = Metric(Array{Float64,4}(params.size), 
                       Array{Float64,4}(params.size), 
                       Array{Float64,4}(params.size), 
                       Array{Float64,4}(params.size), 
                       Array{Float64,4}(params.size), 
                       Array{Float64,4}(params.size))
    for index in CartesianRange(mesh.n)
        SoAmetric.g01[index] = AoSmetric[index].g01 
        SoAmetric.g10[index] = AoSmetric[index].g10 
        SoAmetric.g22[index] = AoSmetric[index].g22 
        SoAmetric.g33[index] = AoSmetric[index].g33 
        SoAmetric.detg[index] = AoSmetric[index].detg 
        SoAmetric.riccis[index] = AoSmetric[index].riccis 
    end
    return SoAmetric
end

function power(x,a)
    return x^a
end

function find_TR_from_UV(U, V)
    @vars x
    u = find_zero(atan((x/(sqrt(1+x^2)))*log(1+x^2)) - U, atan(U))
    v = find_zero(atan((x/(sqrt(1+x^2)))*log(1+x^2)) - V, atan(V))
    t = (u+v)/2
    r = (v-u)/2
    return (t, r)
end

function find_UV_from_TR(t, r)
    u = t-r
    v = t+r
    U = atan((u/(sqrt(1+u^2)))*log(1+u^2))
    V = atan((v/(sqrt(1+v^2)))*log(1+v^2))
    return (U,V)
end

function createmesh(params::Params)::Array{Tuple,2}
    UV = Array{Tuple,2}(params.size)
    for index in CartesianRange(params.size)
        UV[index] = (((params.ulim[2] - params.ulim[1])/2)*chebx(index[1], params.size[1]-1) + sum(params.ulim)/2,
                     ((params.vlim[2] - params.vlim[1])/2)*chebx(index[2], params.size[2]-1) + sum(params.vlim)/2)
    end
    return UV
end
