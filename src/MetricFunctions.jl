#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2018
#--------------------------------------------------------------------

import Base: getindex, size, setindex!, eltype

struct u end
struct d end
struct uu end
struct dd end
struct udd end

struct Metric{Tag, D}
    components::Array{T, 1} where {T}
end

struct Derivative{Tag, D}
    components::Array{T, 1} where {T}
end

mutable struct Christoffel{Tag, D}
    components::Array{Field, 3}
end

mutable struct Ricci{Tag, D}
    components::Array{Field, 1}
end

eltype(r::Metric) = eltype(r.components)
eltype(gamma::Christoffel) = eltype(gamma.components)
eltype(gamma::Ricci) = eltype(gamma.components)

dim(::Metric{Tag, D}) where {Tag, D} = D 

size(g::Metric{Tag, D}, ::Int) where {Tag, D} = D
size(d::Derivative{Tag, D}, ::Int) where {Tag, D} = D
size(d::Christoffel{Tag, D}, ::Int) where {Tag, D} = D
size(d::Ricci{Tag, D}, ::Int) where {Tag, D} = D

Christoffel(g::Metric) = Christoffel{udd, dim(g)}(fill(zero(g[1,1].space), (dim(g), dim(g), dim(g))))
Ricci(g::Metric) = Ricci{dd, dim(g)}(fill(zero(g[1,1].space), (Int(dim(g)*((dim(g)+1)/2)))))

function mapmetricindex(i::Int, j::Int, D::Int)
    return Int(i + D*(j-1) + ((j/2)*(1-j)))
end

function getindex(g::Metric, a::Int, b::Int) 
    #a < b ? (a,b) = (b,a) : (a,b) = (a,b)
    if a >= b
        (i,j) = (a,b)
    else
        (i,j) = (b,a)
    end
    @assert i >= j 
    return g.components[mapmetricindex(i, j, dim(g))] 
end

function getindex(d::Derivative, a::Int) 
    return d.components[a]
end

function getindex(gamma::Christoffel, a::Int, b::Int, c::Int) 
    return gamma.components[a, b, c]
end

function setindex!(gamma::Christoffel, u::Field, a::Int, b::Int, c::Int) 
    gamma.components[a,b,c] = u
    return gamma
end

function getindex(ricci::Ricci, a::Int, b::Int) 
    a < b ? (a,b) = (b,a) : (a,b) = (a,b)
    return ricci.components[mapmetricindex(a, b, size(ricci, 1))]
end

function setindex!(ricci::Ricci, u::Field, a::Int, b::Int) 
    ricci.components[mapmetricindex(a,b, size(ricci, 1))] = u
    return ricci 
end

function metricinverse(g::Metric{dd, 4})::Metric{uu, 4}   
    ginv = Metric{uu, 4}([similar(g[1,1]),
                          similar(g[2,1]), similar(g[2,2]),
                          similar(g[3,1]), similar(g[3,2]), similar(g[3,3]),
                          similar(g[4,1]), similar(g[4,2]), similar(g[4,3]), similar(g[4,4])])

    for index in CartesianRange(size(g[1,1].space))
        tempg = [g[1,1].value[index]  g[1,2].value[index] g[1,3].value[index] g[1,4].value[index];
                 g[2,1].value[index]  g[2,2].value[index] g[2,3].value[index] g[2,4].value[index];
                 g[3,1].value[index]  g[3,2].value[index] g[3,3].value[index] g[3,4].value[index];
                 g[4,1].value[index]  g[4,2].value[index] g[4,3].value[index] g[4,4].value[index]]

        tempginv = inv(tempg)    

        ginv[1,1].value[index] = tempginv[1,1] 
        ginv[2,1].value[index] = tempginv[2,1] 
        ginv[3,1].value[index] = tempginv[3,1] 
        ginv[4,1].value[index] = tempginv[4,1] 
                                             
        ginv[2,2].value[index] = tempginv[2,2] 
        ginv[3,2].value[index] = tempginv[3,2] 
        ginv[4,2].value[index] = tempginv[4,2] 
                                             
        ginv[3,3].value[index] = tempginv[3,3] 
        ginv[4,3].value[index] = tempginv[4,3] 
                                             
        ginv[4,4].value[index] = tempginv[4,4] 

    end
    return ginv
end

function metricdet(g::Metric{dd, 4})::Field 
    detg = similar(g[1,1]) 
    for index in CartesianRange(size(g[1,1].space))
        tempg = [g[1,1].value[index]  g[1,2].value[index] g[1,2].value[index] g[1,4].value[index];
                 g[2,1].value[index]  g[2,2].value[index] g[2,2].value[index] g[2,4].value[index];
                 g[3,1].value[index]  g[3,2].value[index] g[3,2].value[index] g[3,4].value[index];
                 g[4,1].value[index]  g[4,2].value[index] g[4,2].value[index] g[4,4].value[index]]
        detg.value[index] = det(tempg) 
    end
    return detg
end
