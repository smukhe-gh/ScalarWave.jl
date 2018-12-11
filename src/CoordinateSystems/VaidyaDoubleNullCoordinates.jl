#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 06-2018
# Double null coordinates for a Vaidya Spacetime
# See Waugh, Lake 1986
# Roadmap: Choose the definition of m
#          Compute B (eq. 18)
#          Use B to compute r (eq. 12)
#          Compute f from r (eq. 11)
# Test with known analytic solutioin for the linear mass function
#--------------------------------------------------------------------

using ForwardDiff, NLSolve

function computeFR(S::Space, m::Function)
    DV = derivative(S)
    m  = Fiedl(S, u->m(u))
    
    B  = - (DV*m)/2*(abs(DV*m)) 
    A  = m
    # solve the non-linear system 

end

function find_t_of_UV(U::T, V::T, M::Float64)::T where {T<:Float64}
    @assert V > 0   # ensure you're in region I or II
    if U*V == 0     # r = 2M 
        t = 2M*randn()
    elseif U > 0   # r < 2M
        t = -2M*log(U/V)
    elseif U < 0   # r > 2M
        t = -2M*log(-U/V)
    else
        error("Domain error")
    end
    return t
end

function find_r_of_UV(U::T, V::T, M::Float64)::T where {T<:Float64}
    @assert V > 0   # ensure you're in region I or II
    if U*V == 0     # r = 2M 
        r = 2M
    else            # r < 2M or r > 2M 
        f(r) = (r/2M - 1)*exp(r/2M) + U*V 
        r    = find_zero(f, 2M)
    end
    @assert r > 0
    @assert r > 2M      # XXX: For testing.
    return r
end

function find_U_of_tr(t::T, r::T, M::Float64)::T where {T<:Float64}
    @assert r > 0
    if r == 2M
        U = 0
    else
        rstar = r + 2M*log((r/2M)-1)
        u = t - rstar
        r > 2M ?  U = -exp(-u/4M) : U = +exp(-u/4M)
    end
    return U 
end

function find_V_of_tr(t::T, r::T, M::Float64)::T where {T<:Float64}
    @assert r > 0
    if r == 2M
        V = 2M*rand()
    else
        rstar = r + 2M*log((r/2M)-1)
        v = t + rstar
        V = exp(v/4M)
    end
    return V
end


