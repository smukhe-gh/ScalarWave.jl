#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 07-2019
# Test non-linear solver
#--------------------------------------------------------------------

struct U end
struct V end

SU  = ChebyshevGL{U, 2, Float64}(-1, 1)
SV  = ChebyshevGL{V, 2, Float64}( 2, 4)
SUV = ProductSpace(SU, SV)
DU, DV = derivative(SUV)

function stack(u::Field{S}, v::Field{S}, w::Field{S}) where {S}
    return vcat(reshape(u), reshape(v), reshape(w)) 
end

function unstack(space::S, x::Array{T,1}) where {S, T}
    foldx = reshape(x, (:, 3))
    return (Field(space, foldx[:, 1]), 
            Field(space, foldx[:, 2]), 
            Field(space, foldx[:, 3])) 
end

function Base. log(u::Field{S})::Field{S} where {S}
    return Field(u.space, log.(u.value))
end

function f!(F, x)
    (r, f, ϕ) = unstack(SUV, x)
    resr = DU*DV*r + (1/r)*(DU*r)*(DV*r) - f/r # Check this expression since the residual is non-zero
    resf = DU*DV*log(abs(f)) + (2/r)*(DU*DV*r) + 2*(DU*ϕ)*(DV*ϕ)
    resϕ = DU*DV*ϕ + (1/r)*(DU*r)*(DV*ϕ) + (1/r)*(DV*r)*(DU*ϕ)
    F[:] = stack(resr, resf, resϕ)
end

