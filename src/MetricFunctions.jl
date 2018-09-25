#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2018
# Compute the Christoffel symbols and the covariant derivatives
# WORK IN PROGRESS
#--------------------------------------------------------------------

struct metric{Tag, S} 
    guu::Field{S}
    guv::Field{S}
    gvv::Field{S}
end

function christoffels(guu::Field{ProductSpace{S1, S2}}, 
                      guv::Field{ProductSpace{S1, S2}}, 
                      gvv::Field{ProductSpace{S1, S2}}) where {S1, S2}

    g = [guu guv; guv gvv]
    Dv, Du = derivative{ProductSpace{S1, S2}}


end
