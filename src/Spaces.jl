#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2018
# Define operations and data structures for 1D spaces
#--------------------------------------------------------------------

import Base: *, /, +, -, zeros, ones, range, identity

#--- AbstractTypes.jl -----------------------------------------------
abstract type Manifold end
abstract type Space <: Manifold end
abstract type Cardinal <: Space end
#--------------------------------------------------------------------

#--- Spaces.jl ------------------------------------------------------
struct Chebyshev <: Space 
    order::Int
end

struct GaussLobatto <: Cardinal
    order::Int
end

function range(S::Space)
    return 1:S.order+1
end

function dim(S::Space)
    return length(S.order)
end

struct Field{S<:Space, T<:AbstractArray} 
    space::S
    value::T
end

struct Operator{S<:Space, T<:Real}
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
struct ProductSpace <: Space 
    spaces::NTuple{} 
end

struct ProductSpaceOperator <: Space
    space::S where S<:ProductSpace
    value::NTuple{}
end

function dim(S::ProductSpace)
    return length(S.spaces)
end

function ⦼{S<:Cardinal}(S1::S, S2::S)::ProductSpace 
    return ProductSpace((S1, S2))    
end

function ⦼(S1::ProductSpace, S2::Cardinal)::ProductSpace 
    return ProductSpace((S1.spaces..., S2))    
end

function ⦼{S<:Cardinal}(S1::S, S2::S)::ProductSpace 
    return ProductSpace((S1, S2))    
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

function ⦼{A<:Operator}(A1::A, A2::A...)::ProductSpaceOperator 
    return ProductSpaceOperator(ProductSpace((A1.space, A2.space...)), 
                                kron(A1.value, A2.value...))    
end

function derivative(S<:ProductSpace)::ProductSpaceOperator
   orders = map(x->x.order, S.spaces) 
   operators = hcat(derivative, repeat(identity,dim(S)-1))
   D = []
   for combination in unique(collect(permutations((operators), dim(S))),1)
       D = push!(D, ⊙(tuple(map((x, umap)->umap(x), orders, combination)...)...))
   end
   return ProductSpaceOperator(S, tuple(D...))
end

function derivative(S<:ProductSpace)::ProductSpaceOperator
   orders = map(x->x.order, S.spaces) 
   W = ⊙(tuple(map((x, umap)->umap(x), orders, (integral)))...)
   return ProductSpaceOperator(S, (W))
end
#--------------------------------------------------------------------



