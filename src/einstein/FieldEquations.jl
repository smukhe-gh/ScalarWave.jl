#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2019
# See Gundlach and Pullin 1997 for the 
# relevant equations
#--------------------------------------------------------------------

function H(f::Field{S}, r::Field{S}, ϕ::Field{S}, eq::Symbol)::Field{S} where {S}
    DV, DU = derivative(S)
    if eq == :H1
        return (DV*DU)*r + (1/r)*(DU*r)*(DV*r) - (1/r)*f
    elseif eq == :H2
        return DV*DU*log(f) + (2/r)*DV*DU*r + 2*(DU*ϕ)*(DV*ϕ)
    elseif eq == :H3
        return DU*DV*ϕ + ((DU*r)/r)*ϕ + ((DV*r)/r)*ϕ
    else
        @warn "Invalid equation."
        return zero(S) 
    end
end

function E(f::Field{S}, r::Field{S}, ϕ::Field{S}, eq::Symbol)::Field{S} where {S} 
    DV, DU = derivative(S)
    if eq==:E1 
        return DU*DU*r - ((DU*f)/f)*DU*r + (DU*ϕ)^2
    elseif eq==:E2
        return DV*DV*r - ((DV*f)/f)*DV*r + (DV*ϕ)^2
    else
        @warn "Invalid variable"
        return zero(S)
    end
end

function linearH(f::Field{S}, r::Field{S}, ϕ::Field{S}, 
                 eq::Symbol, var::Symbol)::ProductSpaceOperator{S} where {S}
    DV, DU = derivative(S)
    I      = eye(S)
    if eq == :H1
        if var == :Δf
            return -(1/r)*I
        elseif var == :Δr   
            return DU*DV + (1/r)*(DV*r)*DU + (1/r)*(DU*r)*DV + (1/r^2)*f*I - (1/r^2)*(DU*r)*(DV*r)*I
        elseif var == :Δϕ
            return (I - I) 
        else
            @warn "Invalid linearized variable"
        end
    elseif eq == :H2
        if var == :Δf
            return (1/f)*DU*DV - (1/f^2)*(DV*f)*DU - (1/f^2)*(DU*f)*DV + (2/f^3)*(DU*f)*(DV*f)*I - (1/f^2)*(DU*DV*f)*I
        elseif var == :Δr   
            return (2/r)*DU*DV - (1/r^2)*(DU*DV*r)*I
        elseif var == :Δϕ
            return 2*(DU*ϕ)*DV + 2*(DV*ϕ)*DU 
        else
            @warn "Invalid linearized variable"
        end
    elseif eq == :H3
        if var == :Δf
            return (I - I)
        elseif var == :Δr
            return (1/r)*(DV*ϕ)*DU + (1/r)*(DU*ϕ)*DV - (1/r^2)*(DV*r)*(DU*ϕ)*I - (1/r^2)*(DU*r)*(DV*ϕ)*I
        elseif var == :Δϕ
            return DU*DV + (1/r)*(DV*r)*DU + (1/r)*(DU*r)*DV
        else
            @warn "Invalid linearized variable"
        end
    else
        @warn "Invalid equation"
    end
end

function linearE(f::Field{S}, r::Field{S}, ϕ::Field{S}, 
                 eq::Symbol, var::Symbol)::ProductSpaceOperator{S} where {S}
    DV, DU = derivative(S)
    EI     = eye(S)
    if  eq == :E1
        if var == :Δf
            return -(1/f)*(DU*r)*DU + (1/f^2)*(DU*f)(DU*r)*I
        elseif var == :Δr   
            return DU*DU - (1/f)*(DU*f)*DU
        elseif var == :Δϕ
            return 2*(DU*ϕ)*DU
        else
            @warn "Invalid linearized variable"
        end
    elseif eq == :E2
        if var == :Δf
            return -(1/f)*(DV*r)*DV + (1/f^2)*(DV*f)(DV*r)*I
        elseif var == :Δr   
            return DV*DV - (1/f)*(DV*f)*DV
        elseif var == :Δϕ
            return 2*(DV*ϕ)*DV
        else
            @warn "Invalid linearized variable"
        end
    else
        @warn "Invalid equation"
    end
end

function rhsH(f::Field{S}, r::Field{S}, ϕ::Field{S}, eq::Symbol)::Field{S} where {S}
    return -H(f, r, ϕ, eq)
end

function rhsE(f::Field{S}, r::Field{S}, ϕ::Field{S}, eq::Symbol)::Field{S} where {S}
    return -E(f, r, ϕ, eq)
end

function boundaryOP(f::Field{S}, r::Field{S}, ϕ::Field{S}, var::Symbol) where {S}
    B = boundary(S)
    if var == :Δf || var == :Δr || var == :Δϕ
        return B
    else
        return B - B
    end
end

function setBCs!(f::Field{S},   r::Field{S},   ϕ::Field{S},
                 fBC::Field{S}, rBC::Field{S}, ϕBC::Field{S})::Field{S} where {S}
    return ((I-B)*f + B*fBC,
            (I-B)*r + B*rBC,
            (I-B)*ϕ + B*ϕBC)
end
