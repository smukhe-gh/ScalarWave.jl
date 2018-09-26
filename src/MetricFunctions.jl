#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2018
#--------------------------------------------------------------------

import Base: getindex, similar

struct uu end
struct dd end

struct Metric{Tag, D}
    components::Array{T, 1} where {T}
end

dim(::Metric{Tag, D}) where {Tag, D} = D 

function mapmetricindex(i::Int, j::Int, D::Int)
    return Int(i + D*(j-1) + ((j/2)*(1-j)))
end

function getindex(g::Metric, a::Int, b::Int) 
    a < b ? (a,b) = (b,a) : (a,b) = (a,b)
    return g.components[mapmetricindex(a, b, dim(g))]
end

function similar(u::Field{S, D, T})::Field{S,D,T} where {S, D, T}
    return Field(u.space, Array{T,D}(size(u.space)))
end

function metricinverse(g::Metric{dd, 4})::Metric{uu, 4}   
    ginv = Metric{uu, 4}([similar(g[1,1]),
                          similar(g[2,1]), similar(g[2,2]),
                          similar(g[3,1]), similar(g[3,2]), similar(g[3,3]),
                          similar(g[4,1]), similar(g[4,2]), similar(g[4,3]), similar(g[4,4])])

    for index in CartesianRange(size(g[1,1].space))
        tempg = [g[1,1].value[index]  g[1,2].value[index] g[1,2].value[index] g[1,4].value[index];
                 g[2,1].value[index]  g[2,2].value[index] g[2,2].value[index] g[2,4].value[index];
                 g[3,1].value[index]  g[3,2].value[index] g[3,2].value[index] g[3,4].value[index];
                 g[4,1].value[index]  g[4,2].value[index] g[4,2].value[index] g[4,4].value[index]]

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
