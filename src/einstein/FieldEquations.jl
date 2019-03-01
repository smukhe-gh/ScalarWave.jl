#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2019
# See Waugh & Lake 1986 Appendix A for the relevant equations
# TODO: Verify stress-energy tensor terms
# F(x + ΔX) = F(x) + J Δx 
#--------------------------------------------------------------------

function F(f::Field{S}, r::Field{S}, ϕ::Field{S}, eq::Symbol)::Field{S} where {S}
    DV, DU = derivative(S)
    if eq == :θθ
        return (1/f^2)*(r^2)*((1/f)*(DU*f)*(DV*f) - (DU*DV)*f) - (2/f)*r*(DU*DV)*r 
    elseif eq == :UV
        return (1/r^2)*(f + 2*((DU*r)*(DV*r) + r*(DU*DV)*r)) - (DU*ϕ)*(DV*ϕ)
    elseif eq == :UU
        return (2/r)*((1/f)*(DU*r)*(DU*f) - (DU*DU)*r) - (DU*ϕ)^2
    elseif eq == :VV
        return (2/r)*((1/f)*(DV*r)*(DV*f) - (DV*DV)*r) - (DV*ϕ)^2
    elseif eq == :TT
        return DU*DV*ϕ + ((DU*r)/r)*ϕ + ((DV*r)/r)*ϕ
    else
        @warn "Invalid equation."
        return zero(S) 
    end
end

function J( f::Field{S}, r::Field{S}, ϕ::Field{S}, eq::Symbol, 
           Δf::Union{Field{S}, Symbol, Int}, Δr::Union{Field{S}, Symbol, Int}, 
           Δϕ::Union{Field{S}, Symbol, Int})::Union{ProductSpaceOperator{S}, Field{S}} where {S}
    DV, DU = derivative(S)
    if eq == :θθ
        return ( -(1/2)*(r^2/f^2)*(DU*DV*Δf) - (r/f)*(DU*DV*Δr) + (1/2)*(r^2/f^3)*(DV*f)*(DU*Δf) + (1/2)*(r^2/f^3)*(DU*f)*(DV*Δf)
				 -(3/2)*(r^2/f^4)*(DU*f)*(DV*f)*Δf + (r^2/f^3)*(DU*DV*f)*Δf + (r/f^2)*(DU*DV*r)*Δf + (r/f^3)*(DU*f)*(DV*f)*Δr
				 -(r/f^2)*(DU*DV*f)*Δr - (1/f)*(DU*DV*r)*Δr )
    elseif eq == :UV
        return (  (2/r)*DU*DV*Δr + (2/r^2)*(DV*r)*(DU*Δr) + (2/r^2)*(DU*r)*(DV*Δr) + (2/r^2)*Δf 
				 - 4*(f/r^3)*Δr - (4/r^3)*(DU*r)*(DV*r)*Δr - (2/r^2)*(DU*DV*r)*Δr )   
    elseif eq == :UU
        return ( -(2/r)*DU*DU*Δr + (2/f)*(1/r)*(DU*r)*(DU*Δf) + (2/f)*(1/r)*(DU*f)*(DU*Δr) 
				 -(2/f^2)*(1/r)*(DU*f)*(DU*r)*Δf - (2/f)*(1/r^2)*(DU*f)*(DU*r)*Δr
			     +(2/r^2)*(DU*DU*r)*Δr)
    elseif eq == :VV
        return ( -(2/r)*DV*DV*Δr + (2/f)*(1/r)*(DV*r)*(DV*Δf) + (2/f)*(1/r)*(DV*f)*(DV*Δr) 
				 -(2/f^2)*(1/r)*(DV*f)*(DV*r)*Δf - (2/f)*(1/r^2)*(DV*f)*(DV*r)*Δr
			     +(2/r^2)*(DV*DV*r)*Δr )
    elseif eq == :TT
        @error "Stress-energy terms haven't been included yet."
    else
        @warn "Invalid equation"
        return DU*DV - DV*DU
    end
end

function B(space::Type{S}, sym::Symbol)::ProductSpaceOperator{S} where {S} 
    bnd = boundary(Null, space)
    (S == :Δ0) ? (return bnd - bnd) : (return bnd)
end


#--------------------------------------------------------------------
# Construct interfaces to the non-linear solver
#--------------------------------------------------------------------

function Fvec(space::Type{S}, Svec::Array{Float64,1})::Array{Float64,1} where {S}
    (f, r) = Sshape(space, Svec)
    ϕ      = f - f
    return [vec(F(f, r, ϕ, :UV)); 
            vec(F(f, r, ϕ, :θθ))]
end

function Jvec(space::Type{S}, Svec::Array{Float64,1})::Array{Float64,2} where {S}
    (f, r) = Sshape(space, Svec)
    ϕ      = f - f
    return [vec(J(f, r, ϕ, :UV, :Δf, 0, 0)) vec(J(f, r, ϕ, :UV, 0, :Δr, 0)); 
            vec(J(f, r, ϕ, :θθ, :Δf, 0, 0)) vec(J(f, r, ϕ, :θθ, 0, :Δr, 0))] 
end

function Bvec(space::Type{S}, Svec::Array{Float64,1})::Array{Float64,2} where {S}
    return [vec(B(space, :Δf)) vec(B(space, :Δ0)); 
            vec(B(space, :Δ0)) vec(B(space, :Δr))] 
end

function Sshape(space::Type{S}, Svec::Array{Float64,1}) where {S}
    Sshp = reshape(Svec, (prod(size(space)), 2))
    return (Field(space, shape(space, Sshp[:, 1])), 
            Field(space, shape(space, Sshp[:, 2]))) 
end

function Svec(f::Field{S}, r::Field{S}, ϕ::Field{S})::Array{Float64,1} where {S}
    return [vec(f);
            vec(r)]
end
