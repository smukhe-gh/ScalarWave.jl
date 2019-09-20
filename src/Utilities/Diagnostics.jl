#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 04-2019
# Define special operators for axis-symmetry
#--------------------------------------------------------------------

export outgoingexpansion, kretschmannscalar

function outgoingexpansion(a::Field{S}, r::Field{S})::Field{S} where {S}
    DU, DV = derivative(a.space)
    return 2*(DV*r)/(r*abs(-a^2))
end

function kretschmannscalar(a::Field{S}, r::Field{S})::Field{S} where {S}
    # FIXME: Check the expression.
    DU, DV = derivative(a.space)
    f = (1/2)*a^2
    K = (4*(f + 2*(DU*r)*(DV*r))^2/(r^4*f^2) 
         + 16*(f*(DU*(DU*r)) - (DU*r)*(DV*f))*(f*(DV*(DV*r)) - (DV*r)*(DV*f))/(r^2*f^4) 
         + 16*(DU*(DV*r))^2/(r^2*f^2) + 4*((DU*f)*(DV*f) - f*(DU*(DV*f)))^2/(f^6))
    return K
end
