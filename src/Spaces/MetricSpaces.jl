#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2018
# Define metric and derived datatypes along with indexing support.
#--------------------------------------------------------------------

import Base: getindex, size, setindex!, eltype

eltype(g::Metric) = eltype(g.components)
eltype(gamma::Christoffel) = eltype(gamma.components)
eltype(ricci::Ricci) = eltype(ricci.components)
eltype(cd::CovariantDerivative) = eltype(cd.components)

dim(::Metric{Tag, D}) where {Tag, D} = D 

size(g::Metric{Tag, D}, ::Int) where {Tag, D} = D
size(d::Derivative{Tag, D}, ::Int) where {Tag, D} = D
size(d::Christoffel{Tag, D}, ::Int) where {Tag, D} = D
size(d::Ricci{Tag, D}, ::Int) where {Tag, D} = D
size(d::CovariantDerivative{Tag, D}, ::Int) where {Tag, D} = D

Christoffel(g::Metric) = Christoffel{_udd, dim(g)}(fill(zero(g[1,1].space), (dim(g), dim(g), dim(g))))
Ricci(g::Metric) = Ricci{_dd, dim(g)}(fill(zero(g[1,1].space), (Int(dim(g)*((dim(g)+1)/2)))))
CovariantDerivative(g::Metric) = CovariantDerivative{_d, dim(g)}(fill(zero(Null, g.space), Int(dim(g))))

function mapmetricindex(i::Int, j::Int, D::Int)
    return Int(i + D*(j-1) + ((j/2)*(1-j)))
end

function getindex(g::Metric, a::Int, b::Int) 
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

function getindex(CD::CovariantDerivative, a::Int) 
    return CD.components[a]
end

function setindex!(CD::CovariantDerivative, A::ProductSpaceOperator, a::Int) 
    CD.components[a] = A
    return CD 
end
