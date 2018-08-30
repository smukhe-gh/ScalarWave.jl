#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2018
# Define operations and data structures for 1D spaces
#--------------------------------------------------------------------

import Base: size, range, identity 
import Base: +, -, *

struct Chebyshev{Tag, N} <: Galerkin{Tag} end
struct GaussLobatto{Tag ,N} <: Cardinal{Tag} end

struct Field{S, D, T}
    space::Type{S}
    value::Array{T, D}
end

struct Operator{S, D, T}
    space::Type{S}
    value::Array{T, D}
end

order{Tag, N}(S::Type{GaussLobatto{Tag, N}}) = N 
dim{Tag, N}(S::Type{GaussLobatto{Tag, N}})   = 1 
range{Tag, N}(S::Type{GaussLobatto{Tag, N}}) = CartesianRange((N+1,)) 
size{Tag, N}(S::Type{GaussLobatto{Tag, N}})  = range(S).stop.I

+{S}(A::Field{S}, B::Field{S}) = Field(S, A.value + B.value)
-{S}(A::Field{S}, B::Field{S}) = Field(S, A.value - B.value)
*{S}(A::Field{S}, B::Field{S}) = Field(S, A.value .* B.value)
+{S}(A::Operator{S}, B::Operator{S}) = Operator(S, A.value + B.value)
-{S}(A::Operator{S}, B::Operator{S}) = Operator(S, A.value - B.value)

function Field{Tag, N}(S::Type{GaussLobatto{Tag,N}}, umap::Function)::Field{S}
    value = zeros(size(S))
    for index in range(S) 
        value[index] = umap(chebx(index.I[1],N)) 
    end
    return Field(S, value)
end

function derivative{Tag, N}(S::Type{GaussLobatto{Tag,N}})::Operator{S}
    DS = zeros(size(S)..., size(S)...) 
    for index in CartesianRange(size(DS))
        i, j = index.I
        DS[index] = chebd(i, j, N) 
    end
    return Operator(S, DS)
end

function identity{Tag, N}(S::Type{GaussLobatto{Tag,N}})::Operator{S}
    I = zeros(size(S)..., size(S)...) 
    for index in CartesianRange(size(I))
        i ,j = index.I
        I[index] = delta(i,j) 
    end
    return Operator(S, I)
end

function *{S}(A::Operator{S}, B::Operator{S})::Operator{S} 
    C = similar(A.value)
    for index in CartesianRange(size(C))
        i, j = index.I
        C[index] = sum(A.value[i,k]*B[k,j] for k in 1:size(S1))
    end
    return Operator(S, C)
end

function boundary{Tag, N}(S::Type{GaussLobatto{Tag,N}})::Operator{S}
    BS = zeros(size(S)) 
    BS[1] = BS[end] = 1
    return Operator(S, diagm(vec(DS)))
end
