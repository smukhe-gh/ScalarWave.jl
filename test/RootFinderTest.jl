#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 07-2019
# Write a non-linear solver from scratch for Einstein's field
# equations
#--------------------------------------------------------------------

using NLsolve
using Roots

function find_r_of_UV(U::T, V::T, M::T)::T where {T<:Number}
    @assert V > 0       # ensure you're in region I or II
    @assert U*V < 1     # ensure you don't hit the singularity
    if U*V == 0     # r = 2M
        r = 2M
    else            # r < 2M or r > 2M
        f(r) = (r/2M - 1)*exp(r/2M) + U*V
        r    = find_zero(f, (2M, 100M); verbose=true)
    end
    @assert r > 0
    return r
end

T = BigFloat
(U, V) = (-4.223542864692437712953303487127855014952793296040208458557990378054829075753572, 6.73205080756887729352744634150587236694280525381038062805580697945193301690878)
U = T(U)
V = T(V) 
M = T(1)

@show find_r_of_UV(U, V, M)
