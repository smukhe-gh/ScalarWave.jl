#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
# Rationals Test
# Choose uniform collocation points in a cardinal basis (polynomial)
# and compute the derivative and the integral
#--------------------------------------------------------------------

function collocation(::Type{Rational}, i::Int, N::Int)::Rational
    return -(-1 + (2*(i-1)//N))
end

function collocation(::Type{Float64}, i::Int, N::Int)::Float64
    return chebx(i, N) 
end

function deriv(::Type{Rational}, i::Int, j::Int, N::Int)::Rational
    x = [collocation(Rational, k, N) for k in 1:N+1]
    if i == j
        return sum((k==j ? 0 : 1/(x[j] - x[k])) for k in 1:N+1)
    else
        ai = prod((k==i ? 1 : (x[i] - x[k])) for k in 1:N+1)
        aj = prod((k==j ? 1 : (x[j] - x[k])) for k in 1:N+1)
        return ai/(aj*(x[i] - x[j]))
    end
end

function deriv(::Type{Float64}, i::Int, j::Int, N::Int)::Float64
    return chebd(i, j, N)
end

