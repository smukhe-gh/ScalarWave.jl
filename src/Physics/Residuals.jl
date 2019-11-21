#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 11-2019
# Solve constraint equations on in the initial hypersurface
# We use a variable rescaling similar to Garfinkle.
#--------------------------------------------------------------------

export C1, C2, E1, E2, E3

function C1(s::Field{S}, ψ::Field{S}, ϕ::Field{S})::Field{S} where {S<:ProductSpace{S1, S2}} where {S1, S2}
    DU, DV = derivative(s.space)
    r = Field(s.space, (u,v)->(v-u))
    return (4*π*(DU*ϕ)^2 - 4*(r^2)*ψ*(DU*s)*(DU*ψ) 
            - 2*s*ψ*(DU*r)*(ψ*(DU*r) + 2*r*(DU*ψ)) 
            + (ψ^2)*(DU*(DU*r)) - 2*r*(ψ^2*(DU*r)*(DU*s) + 3*(DU*ψ)^2 - ψ*(DU*(DU*ψ))))
end

function C2(s::Field{S}, ψ::Field{S}, ϕ::Field{S})::Field{S} where {S<:ProductSpace{S1, S2}} where {S1, S2}
    DU, DV = derivative(s.space)
    r = Field(s.space, (u,v)->(v-u))
    return (4*π*(DV*ϕ)^2 - 4*(r^2)*ψ*(DV*s)*(DV*ψ) 
            - 2*s*ψ*(DV*r)*(ψ*(DV*r) + 2*r*(DV*ψ)) 
            + (ψ^2)*(DV*(DV*r)) - 2*r*(ψ^2*(DV*r)*(DV*s) + 3*(DV*ψ)^2 - ψ*(DV*(DV*ψ))))
end

function E1(s::Field{S}, ψ::Field{S}, ϕ::Field{S})::Field{S} where {S<:ProductSpace{S1, S2}} where {S1, S2}
    DU, DV = derivative(s.space)
    r = Field(s.space, (u,v)->(v-u))
    return ((1/r)*(ψ^2)*exp(2*r*s) + (1/r)*(ψ^2)*(DU*r)*(DV*r) 
            + 4*ψ*(DV*r)*(DU*ψ) + 4*ψ*(DU*r)*(DV*ψ) 
            + 6*r*(DU*ψ)*(DV*ψ) + (ψ^2)*(DU*(DV*r)) + 2*r*ψ*(DU*(DV*ψ)))
end

function E2(s::Field{S}, ψ::Field{S}, ϕ::Field{S})::Field{S} where {S<:ProductSpace{S1, S2}} where {S1, S2}
    DU, DV = derivative(s.space)
    r = Field(s.space, (u,v)->(v-u))
    return (2*(DV*r)*(DU*s) + 2*(DU*r)*(DV*s)
            + 4*π*(DU*ϕ)*(DV*ϕ) + (4/(r*ψ))*((DV*r)*(DU*ψ) + (DU*r)*(DV*ψ))
            + (2/r)*(DU*(DV*r)) + 2*r*(DU*(DV*s)) + (8/ψ)*(DU*(DV*ψ)))
end

function E3(s::Field{S}, ψ::Field{S}, ϕ::Field{S})::Field{S} where {S<:ProductSpace{S1, S2}} where {S1, S2}
    DU, DV = derivative(s.space)
    r = Field(s.space, (u,v)->(v-u))
    return ((1/r)*(DV*r)*(DU*ϕ) + (1/r)*(DU*r)*(DV*ϕ) 
            + (2/ψ)*(DV*ϕ)*(DU*ψ) + (2/ψ)*(DU*ϕ)*(DV*ψ) + (DU*(DV*ϕ)))
end
