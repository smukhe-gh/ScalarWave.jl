#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2018
# Add basis functions
#--------------------------------------------------------------------

export Cardinal, PointSpace
export ChebyshevGL, LegendreGL, FourierEP, Chebyshev

struct ChebyshevGL{Tag, N, T} <: Space{Tag, N} 
    min::T
    max::T
    ChebyshevGL{Tag, N, T}(min, max) where {Tag, N, T} = max > min ? new{Tag, N, T}(min, max) : error("bounds are out of order")
end

struct LegendreGL{Tag, N, T} <: Space{Tag, N} 
    min::T
    max::T
    LegendreGL{Tag, N, T}(min, max) where {Tag, N, T} = max > min ? new{Tag, N, T}(min, max) : error("bounds are out of order")
end

struct FourierEP{Tag, N, T} <: Space{Tag, N} 
    min::T
    max::T
    FourierEP{Tag, N, T}(min, max) where {Tag, N, T} = new{Tag, N, T}(0, 2*pi) 
end

Cardinal{Tag, N, T} = Union{ChebyshevGL{Tag, N, T}, LegendreGL{Tag, N, T}, FourierEP{Tag, N, T}} 

struct Chebyshev{Tag, N, T} <: Space{Tag, N} 
    min::T
    max::T
    Chebyshev{Tag, N, T}(min, max) where {Tag, N, T} = max > min ? new{Tag, N, T}(min, max) : error("bounds are out of order")
end

