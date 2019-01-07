#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2019
#--------------------------------------------------------------------

using Calculus, Einsum

struct SymField 
    expression::Expr
end

struct SymDerivative
    parameter::Symbol
end

struct ProductSymDerivative
    parameters::Array{Symbol}
end

struct MixedDerivative
    parameters::Array{Union{SymDerivative, ProductSpaceOperator}}
end

function Base. *(D::SymDerivative, E::SymDerivative)
    return ProductSymDerivative([D.parameter, E.parameter])
end

function Base. *(D::ProductSpaceOperator, E::SymDerivative)
    return MixedDerivative([E, D])
end

function Base. *(D::SymDerivative, E::ProductSpaceOperator)
    return MixedDerivative([D, E])
end

function Base. *(D::SymDerivative, u::SymField)
    expr = simplify(differentiate(u.expression, D.parameter))
    if expr == 0
        return 0
    else
        return SymField(simplify(differentiate(u.expression, D.parameter)))
    end
end

function Base. *(D::ProductSymDerivative, u::SymField)
    expr = u.expression
    for symbol in D.parameters
        expr = simplify(differentiate(expr, symbol))
    end
    if expr == 0
        return 0
    else
        return SymField(simplify(expr))
    end
end

function Base. *(D::MixedDerivative, u::SymField)
    expr = u
    for operator in D.parameters
        expr = operator*expr 
    end
    return expr
end

function Base. *(D::ProductSpaceOperator, u::SymField)
    return D*eval(u.expression)
end

function Base. *(D::T, u::Field) where {T<:Union{ProductSymDerivative, SymDerivative, MixedDerivative}}
    @warn "A symbolic derivative on a field returns zero by default.
           Use this operation with caution."
    return 0
end

function Base. +(u::Field, a::Int64)
    return a + u
end

#-------------------------------------------------------------------
# Test case
#-------------------------------------------------------------------

SVU = ProductSpace{GaussLobatto(V, 3, 4, 2), 
                   GaussLobatto(U, 3, 4, 2)} 

Dr, Dt = ScalarWave.derivative(SVU)
Dθ, Dϕ = SymDerivative(:θ), SymDerivative(:ϕ)

M = 1
r = Field(SVU, (u,v) -> u+v)
θ = Field(SVU, (u,v) -> u*v)

gtt = (1-2M/r)^(-1) 
grr = (1-2M/r)
gθθ = r^2
gϕϕ = SymField(:(r^2*sin(θ)^2))

function Base. +(u::Field, v::SymField)
    return u + eval(v.expression)
end

# compute the Christoffels
# and then the Riemann curvature tensor

g = [gtt  0   0   0; 
      0  grr  0   0; 
      0   0  gθθ  0; 
      0   0   0  gϕϕ]

D = [Dt, Dr, Dθ, Dϕ]

#@show sum(D[m]*g[m,m] for m in 1:4)

@show Dθ*Dθ*gϕϕ
@show Dr*Dθ*gϕϕ
@show Dθ*Dr*gϕϕ
@show Dr*Dr*gϕϕ

@show Dθ*Dθ*grr
@show Dr*Dθ*grr
@show Dθ*Dr*grr
@show Dr*Dr*grr
