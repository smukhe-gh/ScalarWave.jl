#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Coordinate transformatins for (compactified) 
# Schwarzschild background in double-null coordinates
# Soham 06-2018
#--------------------------------------------------------------------

function find_TR_of_UV(U::Float64, V::Float64, M::Float64)::Tuple
    @assert U*V < 1
    f(x) = exp(x/2M)*(x/2M - 1) + U*V
    r = find_zero(f, (2*M))
    t = -2M*log(abs(U/V))
    return (t, r)
end

function find_UV_of_TR(t::Float64, r::Float64, M::Float64)::Tuple
    @assert r > 0
    @assert r > 2M
    rs = r + 2M*log(r/2M-1)
    u  = t-rs
    v  = t+rs
    U  = r > 2M ? -exp(-u/4M) : exp(-u/4M) 
    V  = exp(v/4M)
    @assert U*V < 1
    return (U,V)
end
