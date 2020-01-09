#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2019
# Compute Initial Data for rescaled variables 
#--------------------------------------------------------------------

using NLsolve, ForwardDiff
export L_for_s, rhs_for_s, s_on_ubnd, rescale, unscale, rescaledC2onbnd, nonlinearsolver_for_psi

function L_for_s(ψ::Field{S}, ϕ::Field{S}, r::Field{S}, DV::Operator{S}, I::Operator{S})::Operator{S} where {S<:Space{Tag}} where {Tag}
    return -2*(ψ^2)*((DV*r)^2)*I - 2*r*(ψ^2)*(DV*r)*DV - 4*r*ψ*(DV*r)*(DV*ψ)*I - 4*(r^2)*ψ*(DV*ψ)*DV
end

function rhs_for_s(ψ::Field{S}, ϕ::Field{S}, r::Field{S}, DV::Operator{S})::Field{S} where {S<:Space{Tag}} where {Tag}
    return (4*π*(DV*ϕ)^2 - 6*r*(DV*ψ)^2 + (ψ^2)*(DV*(DV*r)) + 2*r*ψ*(DV*(DV*ψ)))
end

function s_on_ubnd(s::Field{S}, ψ::Field{S}, ϕ::Field{S}) where {S<:ProductSpace{S1, S2}} where {S1, S2}
    r = Field(s.space, (u,v)->(v-u))
    ronUbnd = extractUboundary(r, :incoming)
    sonUbnd = extractUboundary(s, :incoming)
    ψonUbnd = extractUboundary(ψ, :incoming)
    ϕonUbnd = extractUboundary(ϕ, :incoming)
    DV = derivative(ronUbnd.space)
    I = identity(ronUbnd.space)
    B = incomingboundary(ronUbnd.space)

    L = L_for_s(ψonUbnd, ϕonUbnd, ronUbnd, DV, I)
    b = rhs_for_s(ψonUbnd, ϕonUbnd, ronUbnd, DV)
    
    @show cond(L ⊕ B)
    sonUbnd = solve(L ⊕ B, (B*sonUbnd) ⊕ (-b))
    @show L2(L*sonUbnd + b)
    @show ronUbnd
    @show (L*sonUbnd + b)

    return sonUbnd
end

function rescale(a::Field{S}, η::Field{S}, ϕ::Field{S})::NTuple{3, Field{S}} where {S}
    r = Field(a.space, (u,v)->(v-u))
    f = a^2
    ψ = sqrt(η/r)
    s = log(f/(2*ψ^4))/(2*r)
    return (s, ψ, ϕ)
end

function unscale(s::Field{S}, ψ::Field{S}, ϕ::Field{S})::NTuple{3, Field{S}} where {S}
    r = Field(s.space, (u,v)->(v-u))
    f = 2*(ψ^4)*exp(2*r*s)
    a = sqrt(f)
    η = r*ψ^2 
    return (a, η, ϕ)
end

function rescaledC2onbnd(s::Field{S}, ψ::Field{S}, ϕ::Field{S}, r::Field{S}, DV::Operator{S})::Field{S} where {S<:Space{Tag}} where {Tag}
    return (4*π*(DV*ϕ)^2 - 4*(r^2)*ψ*(DV*s)*(DV*ψ) 
        - 2*s*ψ*(DV*r)*(ψ*(DV*r) + 2*r*(DV*ψ)) 
        + (ψ^2)*(DV*(DV*r)) - 2*r*(ψ^2*(DV*r)*(DV*s) + 3*(DV*ψ)^2 - ψ*(DV*(DV*ψ))))
end

function nonlinearsolver_for_psi(s::Field{S}, ψ::Field{S}, ϕ::Field{S})::Field{S2} where {S<:ProductSpace{S1, S2}} where {S1, S2}
    r = Field(s.space, (u,v)->(v-u))
    
    # Convert to data on boundary
    r = extractUboundary(r, :incoming)
    s = extractUboundary(s, :incoming)
    ψ = extractUboundary(ψ, :incoming)
    ϕ = extractUboundary(ϕ, :incoming)

    DV = derivative(r.space)
    I = identity(r.space)
    B  = incomingboundary(r.space) + outgoingboundary(r.space)
    ψbnd = B*ψ
    
    function F(ψ::Field{S})::Field{S} where {S}
        res = (4*π*(DV*ϕ)^2 - 4*(r^2)*ψ*(DV*s)*(DV*ψ) 
           - 2*s*ψ*(DV*r)*(ψ*(DV*r) + 2*r*(DV*ψ)) 
           + (ψ^2)*(DV*(DV*r)) - 2*r*(ψ^2*(DV*r)*(DV*s) + 3*(DV*ψ)^2 - ψ*(DV*(DV*ψ))))
        return (I - B)*res + B*(ψ - ψbnd)
    end
    
    function f!(fvec::Array{T,1}, x::Array{T,1}) where {T}
        fvec[:] = reshape(F(reshape(r.space, x)))
    end
    
    function j!(jvec::Array{T,2}, x::Array{T,1}) where {T}
        jvec[:, :] = ForwardDiff.jacobian(f!, similar(x), x)
    end
    
    ψsolved = reshape(r.space, nlsolve(f!, reshape(r + ψ); method=:trust_region, autodiff=:forward, show_trace=true, ftol=1e-10, iterations=120).zero)
    return ψsolved
end
