#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2018
# Types for 1D setup [Add dimensions whereever necessary]
#--------------------------------------------------------------------

struct Manifold{T<:Real}
    min::T
    max::T
    npoints::Int
end

struct RealSpace{T<:Real}
    manifold::Manifold
    class::String
    elements::Array{T,1}
end

function chart(M::Manifold, class::String)::RealSpace 
    A = (M.max + M.min)/2
    B = (M.max - M.min)/2
    points = map(x->B*chebx(x, M.npoints-1) + A, range(1, M.npoints))
    return RealSpace(M, class, points)    
end

struct Field{T<:Real}
    manifold::Manifold
    elememts::Array{T,1}
end

function *(u::Field, v::Field)::Field
    @assert u.manifold == v.manifold
    return Field(u.manifold, u.elements.*v.elements) 
end

function +(u::Field, v::Field)::Field
    @assert u.manifold == v.manifold
    return Field(u.manifold, u.elements.+v.elements) 
end

struct VectorSpace{T<:Real}
    # NOTE: How do we define dimension?
    field::Field 
    dimension::Int
    elements::Array{T,1}
end

struct TensorSpace{T<:Real}
    u::VectorSpace
    v::VectorSpace
    dimension::Int
    elements::Array{T,2}
end

function +(u::VectorSpace, v::VectorSpace)::VectorSpace
    # NOTE: Can two vectorspaces defined on two different fields 
    #       be added?
    @assert u.dimension == v.dimension
    return VectorSpace(u.field, u.dimension, u.elements .+ v.elements)
end

function *(u::Real, v::VectorSpace)::VectorSpace
    # NOTE: Can two vectorspaces defined on two different fields 
    #       be added?
    @assert u.dimension == v.dimension
    return VectorSpace(u.field, u.dimension, u.*v.elements)
end

function ⦼(u::VectorSpace, v::VectorSpace)::TensorSpace
    # FIXME: Needs to be generalized to higher dimensions
    #        Can vectorspaces be clubbed into tensor spaces?
    return TensorSpace(u, v, u.dimension*v.dimension, 
                            kron(u.elements, v.elements))
end

struct DualSpace{T<:Real}
    dualto::VectorSpace
    elements::Array{T,1}
end

function *(u::DualSpace, v::VectorSpace)::Real
    @assert u.dualto == v
    return sum(u.elements.*v.elements)
end

struct Operator{T<:Real}
    # TODO: We want to define this in terms of co-vectors
    #       living in the dualSpace.
    space::Field
    elements::Array{T}{VectorSpace.dimension, VectorSpace.dimension}
end

function derivOperator(u::Field)::Operator
    # TODO: The derivative operator needs to know about the chart. 
    #       However, there's no natural way to do this, since the
    #       value of the field doens't depend on the chart.
    D = Float64[chebd(i, j, u.manifold.npoints-1) 
                for i in range(1, u.manifold.npoints),
                    j in range(1, u.manifold.npoints)]
    return Operator(u, D)
end

function *(D::Operator, u::Field)::Field
    # FIXME: Generalize this to work with 
    #        arbitrary number of dimensions 
    # Also, this is ugly. 
    return Field(u.manifold, D.elements*u.elements)
end


