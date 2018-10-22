#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2018
# Coordinate transforms for the metric and the derivatives.
#--------------------------------------------------------------------

function inversemetrictransform(guu::T, guv::T, gvv::T, 𝒖::T, 𝒗::T) where {T<:Field{ProductSpace{S1, S2}}}  where {S1, S2}
    𝔻v, 𝔻u = derivative(ProductSpace{S1, S2})
    𝔻uof𝒖 = 𝔻u*𝒖
    𝔻vof𝒖 = 𝔻v*𝒖
    𝔻uof𝒗 = 𝔻u*𝒗
    𝔻vof𝒗 = 𝔻v*𝒗

    DxofX = [𝔻uof𝒖 𝔻vof𝒖; 𝔻uof𝒗 𝔻vof𝒗]
    ginv  = [guu guv; guv gvv]
    𝐠     = [sum(DxofX[a,m]*DxofX[b,n]*ginv[m,n] for m in 1:2, n in 1:2) for a in 1:2, b in 1:2] 

    return (𝐠[1,1], 𝐠[1,2], 𝐠[2,1], 𝐠[2,2])
end

function derivativetransform(PS::Type{ProductSpace{S1, S2}}, 𝒖::Field{ProductSpace{S1, S2}}, 
                                                             𝒗::Field{ProductSpace{S1, S2}}) where {S1, S2}
    𝔻v, 𝔻u = derivative(ProductSpace{S1, S2})

    𝔻uof𝒖 = 𝔻u*𝒖
    𝔻vof𝒖 = 𝔻v*𝒖
    𝔻uof𝒗 = 𝔻u*𝒗
    𝔻vof𝒗 = 𝔻v*𝒗
    
    𝔻𝒖ofu = Field(ProductSpace{S1, S2}, similar(𝔻uof𝒖.value)) 
    𝔻𝒖ofv = Field(ProductSpace{S1, S2}, similar(𝔻vof𝒖.value)) 
    𝔻𝒗ofu = Field(ProductSpace{S1, S2}, similar(𝔻uof𝒗.value))
    𝔻𝒗ofv = Field(ProductSpace{S1, S2}, similar(𝔻vof𝒗.value))
    
    for index in CartesianRange(size(𝔻uof𝒖.value)) 
        # TODO: This can be made faster, by moving array 
        #       allocation outside the loop.
        Jacobian = [𝔻uof𝒖.value[index] 𝔻uof𝒗.value[index]; 
                    𝔻vof𝒖.value[index] 𝔻vof𝒗.value[index]]
        InverseJacobian    = inv(Jacobian)
        𝔻𝒖ofu.value[index] = InverseJacobian[1,1] 
        𝔻𝒖ofv.value[index] = InverseJacobian[1,2] 
        𝔻𝒗ofu.value[index] = InverseJacobian[2,1] 
        𝔻𝒗ofv.value[index] = InverseJacobian[2,2] 

        DxofX = [𝔻uof𝒖.value[index] 𝔻uof𝒗.value[index];
                 𝔻vof𝒖.value[index] 𝔻vof𝒗.value[index]]

        DXofx = [𝔻𝒖ofu.value[index] 𝔻𝒖ofv.value[index];
                 𝔻𝒗ofu.value[index] 𝔻𝒗ofv.value[index]]
    end

    𝔻𝒖    = 𝔻𝒖ofu * 𝔻u + 𝔻𝒖ofv * 𝔻v  
    𝔻𝒗    = 𝔻𝒗ofu * 𝔻u + 𝔻𝒗ofv * 𝔻v
    
    return(𝔻𝒗, 𝔻𝒖)
end
