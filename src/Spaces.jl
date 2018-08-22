#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2018
# Define operations and data structures for 1D spaces
#--------------------------------------------------------------------

import Base: *, /, +, -, zeros, ones, range, identity

#--- AbstractTypes.jl -----------------------------------------------
abstract type Manifold{Tag} end
abstract type Space{Tag} <: Manifold{Tag} end
abstract type Cardinal{Tag} <: Space{Tag} end
#--------------------------------------------------------------------

#--- Spaces.jl ------------------------------------------------------
struct Chebyshev{Tag} <: Space{Tag} 
    order::Int
end
order(s::Chebyshev) = s.order

struct ChebyshevN{Tag, N} <:Space{Tag} end
order{Tag, N}(::ChebyshevN{Tag, N}) = N

struct GaussLobatto{Tag} <: Cardinal{Tag}
    order::Int
end

function range(S::Space)
    return 1:S.order+1
end

function dim(S::Space)
    return length(S.order)
end

dim(s::Chebyshev) = 1
dim(s::ProductSpace{S1, S2}) = dim(s.space1) + dim(s.space2)
dim(s::UnitSpace) = 0

dim(s::Space) = length(range(s).start)

range(s::Chebyshev) = CartesianRange(1:order(s))
function range(s::ProductSpace(S1, S2))
    start = CartesianIndex(range(s.space1).start, range(s.space2).start)
    stop = CartesianIndex(range(s.space1).stop, range(s.space2).stop)
    CartesianRange(start, stop)
end
range(s::UnitSpace) = CartesianRange(())

function range(S::Space)

struct Field{S<:Space, T<:AbstractArray} 
    space::S
    value::T
end

struct Field{S<:Space, D, T}
    space::S
    value::Array{T, D}
    function Field{S, D, T}(s)
        @assert D == dim(s)
        new{S, D, T}(s, Array{T}(range(s)))
    end
end

struct Operator{S<:Space, T}
    space::S
    value::Array{T}
end

Field(S::Cardinal) = Field(S, chebgrid(S.order)) 
Field(S::Cardinal, umap::Function) = Field(S, umap.(chebgrid(S.order))) 

+{S}(u::Field{S}, v::Field{S})::Field{S} = Field(u.space, u.value .+ v.value) 
-{S}(u::Field{S}, v::Field{S})::Field{S} = Field(u.space, u.value .- v.value) 
*{S<:Cardinal}(u::Field{S}, v::Field{S})::Field{S} = Field(u.space, u.value .* v.value) 
/{S<:Cardinal}(u::Field{S}, v::Field{S})::Field{S} = Field(u.space, u.value ./ v.value) 

zeros{S<:Space}(space::S)::Field{S} = Field(space, zeros(space.order+1))
ones{S<:Space}(space::S)::Field{S}  = Field(space, zeros(space.order+1))

+{S}(A::Operator{S}, B::Operator{S})::Operator{S} = Operator(A.space, A.value .+ B.value)
-{S}(A::Operator{S}, B::Operator{S})::Operator{S} = Operator(A.space, A.value .- B.value)
*{S}(A::Operator{S}, B::Operator{S})::Operator{S} = Operator(A.space, A.value*B.value)
*{S}(A::Operator{S}, u::Field{S})::Field{S} = Field(A.space, A.value*u.value)
*{S<:Cardinal}(u::Field{S}, A::Operator{S})::Operator{S} = Operator(A.space, (eye(u.space.order+1).*u.value)*A.value)

function identity(S::GaussLobatto)::Operator
    return Operator(S, eye(S.order+1))    
end

function derivative(S::GaussLobatto)::Operator
    D = [chebd(i, j, S.order) for i in range(S), j in range(S)]
   return Operator(S, D)
end

function integral(S::GaussLobatto)::Operator
    W = [chebw(i, S.order) for i in range(S)]
    return Operator(S, diagm(W))
end

function boundary(S::GaussLobatto)::Operator
   B = zeros(S.order+1, S.order+1)
   B[1,1]      = 1/chebw(1, S.order)
   B[end, end] = 1/chebw(S.order+1, S.order)
   return Operator(S, B)
end

function solve{S<:Space}(A::Operator{S}, u::Field{S})::Field{S}
    return Field(A.space, A.value\u.value) 
end
#--------------------------------------------------------------------

#--- ProductSpaces.jl -----------------------------------------------
struct ProductSpace{S1<:Space, S2<:Space} <: Space
    space1::S1
    space2::S2
end

struct ProductSpace{S<:NTuple{D,Space}} <: Space 
    spaces::S
end

struct ProductSpace <: Space
    spaces::NTuple{}
end

struct UnitSpace{Tag} <: Space end

struct ProductSpaceOperator <: Space
    space::S where S<:ProductSpace
    value::NTuple{}
end

function dim(S::ProductSpace)
    return length(S.spaces)
end

dim(s::UnitSpace) = 0

function ⦼{S<:Cardinal}(S1::S, S2::S)::ProductSpace 
    return ProductSpace((S1, S2))    
end

function ⦼(S1::ProductSpace, S2::Cardinal)::ProductSpace 
    return ProductSpace((S1.spaces..., S2))    
end

function chebgrid(N::Tuple)::Array{NTuple} 
    X = Array{NTuple}(N.+1)
    for index in CartesianRange(size(X))
        X[index] = map((i, n)->chebx(i, n), index.I, N)
    end
    return X
end

Field(S::ProductSpace) = Field(S, chebgrid(map(s->s.order, S.spaces))) 
Field(S::ProductSpace, umap::Function) = Field(S, map(x->umap.(x...), Field(S).value))

function ⦼(A1::Operator, A2::Operator...)::ProductSpaceOperator 
    spaces =      (A1.space, map(x->x.space, A2)...)
    value  = kron((A1.value, map(x->x.value, A2)...)..)
    return ProductSpaceOperator(ProductSpace(spaces),
                                tuple(value...))
end

function derivative(S::ProductSpace)::ProductSpaceOperator
   order    = map(x->x.order, S.spaces) 
   operator = hcat(derivative, repeat(identity, dim(S)-1))
   D        = Array{ProductSpaceOperator, 1}

   for (index, permutation) in enumerate(unique(collect(permutations((operators), dim(S))),1))
       D[index] = ⊙(tuple(map((S, operator)->operator(S), spaces, permutation)...)...)
   end

   return ProductSpaceOperator(S, tuple(D...))
end

function integral(S::ProductSpace)::ProductSpaceOperator
   orders = map(x->x.order, S.spaces) 
   W = ⊙(tuple(map((S, operator)->operator(S), spaces, [integral])...)...)   
   return ProductSpaceOperator(S, tuple(W...))
end
#--------------------------------------------------------------------



