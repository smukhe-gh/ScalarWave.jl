#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Coordinate transformatins for (compactified) 
# Schwarzschild background in double-null coordinates
# Soham 06-2018
#--------------------------------------------------------------------

function find_t_of_uv(u::T, v::T, M::T)::T where {T<:Real}
    @show u*v
    @assert u*v < 1
    f(x) = exp(x/2M)*(x/2M - 1) + u*v
    t = -2M*log(abs(u/v))
    return t
end

function find_r_of_uv(u::T, v::T, M::T)::T where {T<:Real}
    @assert u*v < 1
    f(x) = exp(x/2M)*(x/2M - 1) + U*V
    r = find_zero(f, (2*M))
    return r
end

function find_u_of_tr(t::T, r::T, M::T)::T where {T<:Real}
    @assert r > 0
    @assert r > 2M
    rs = r + 2M*log(r/2M-1)
    u  = t-rs
    u  = r > 2M ? -exp(-u/4M) : exp(-u/4M) 
    return u 
end

function find_v_of_tr(t::T, r::T, M::T)::T where {T<:Real}
    @assert r > 0
    @assert r > 2M
    rs = r + 2M*log(r/2M-1)
    v  = t+rs
    v  = exp(v/4M)
    return  v
end
