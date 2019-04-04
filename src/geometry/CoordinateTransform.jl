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

function inversemetrictransform(guu::T, guv::T, guw::T, 
                                        gvv::T, gvw::T,
                                                gww::T, 𝒖::T, 𝒗::T, 𝐰::T) where {T<:Field{ProductSpace{S1, S2, S3}}}  where {S1, S2, S3}
    𝔻w, 𝔻v, 𝔻u = derivative(ProductSpace{S1, S2, S3})

    𝔻uof𝒖 = 𝔻u*𝒖
    𝔻vof𝒖 = 𝔻v*𝒖
    𝔻wof𝒖 = 𝔻w*𝒖

    𝔻uof𝒗 = 𝔻u*𝒗
    𝔻vof𝒗 = 𝔻v*𝒗
    𝔻wof𝒗 = 𝔻w*𝒗

    𝔻uof𝐰 = 𝔻u*𝐰
    𝔻vof𝐰 = 𝔻v*𝐰
    𝔻wof𝐰 = 𝔻w*𝐰

    DxofX = [𝔻uof𝒖 𝔻vof𝒖 𝔻wof𝒖; 
             𝔻uof𝒗 𝔻vof𝒗 𝔻wof𝒗;
             𝔻uof𝐰 𝔻vof𝐰 𝔻wof𝐰]

    ginv  = [guu guv guw; 
             guv gvv gvw;
             guw gvw gvv]

    𝐠     = [sum(DxofX[a,m]*DxofX[b,n]*ginv[m,n] for m in 1:3, n in 1:3) for a in 1:3, b in 1:3] 

    return (𝐠[1,1], 𝐠[1,2], 𝐠[2,3], 
            𝐠[1,2], 𝐠[2,2], 𝐠[2,3],
                            𝐠[3,3])
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
    
    for index in CartesianIndices(size(𝔻uof𝒖.value)) 
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

function derivativetransform(PS::Type{ProductSpace{S1, S2, S3}}, 𝒖::Field{ProductSpace{S1, S2, S3}}, 
                                                                 𝒗::Field{ProductSpace{S1, S2, S3}},
                                                                 𝐰::Field{ProductSpace{S1, S2, S3}}) where {S1, S2, S3}
    𝔻v, 𝔻u = derivative(ProductSpace{S1, S2, S3})

    𝔻uof𝒖 = 𝔻u*𝒖
    𝔻vof𝒖 = 𝔻v*𝒖
    𝔻wof𝒖 = 𝔻w*𝒖

    𝔻uof𝒗 = 𝔻u*𝒗
    𝔻vof𝒗 = 𝔻v*𝒗
    𝔻wof𝒗 = 𝔻w*𝒗
    
    𝔻uof𝐰 = 𝔻u*𝐰
    𝔻vof𝐰 = 𝔻v*𝐰
    𝔻wof𝐰 = 𝔻w*𝐰

    𝔻𝒖ofu = Field(ProductSpace{S1, S2, S3}, similar(𝔻uof𝒖.value)) 
    𝔻𝒖ofv = Field(ProductSpace{S1, S2, S3}, similar(𝔻vof𝒖.value)) 
    𝔻𝒖ofw = Field(ProductSpace{S1, S2, S3}, similar(𝔻wof𝒖.value)) 

    𝔻𝒗ofu = Field(ProductSpace{S1, S2, S3}, similar(𝔻uof𝒗.value))
    𝔻𝒗ofv = Field(ProductSpace{S1, S2, S3}, similar(𝔻vof𝒗.value))
    𝔻𝒗ofw = Field(ProductSpace{S1, S2, S3}, similar(𝔻wof𝒗.value))
    
    𝔻𝐰ofu = Field(ProductSpace{S1, S2, S3}, similar(𝔻uof𝐰.value))
    𝔻𝐰ofv = Field(ProductSpace{S1, S2, S3}, similar(𝔻vof𝐰.value))
    𝔻𝐰ofw = Field(ProductSpace{S1, S2, S3}, similar(𝔻wof𝐰.value))

    for index in CartesianIndices(size(𝔻uof𝒖.value)) 
        # TODO: This can be made faster, by moving array 
        #       allocation outside the loop.
        Jacobian = [𝔻uof𝒖.value[index] 𝔻uof𝒗.value[index] 𝔻uof𝐰.value[index];
                    𝔻vof𝒖.value[index] 𝔻vof𝒗.value[index] 𝔻vof𝐰.value[index];
                    𝔻wof𝒖.value[index] 𝔻wof𝒗.value[index] 𝔻wof𝐰.value[index]]
        
        InverseJacobian    = inv(Jacobian)
        𝔻𝒖ofu.value[index] = InverseJacobian[1,1] 
        𝔻𝒖ofv.value[index] = InverseJacobian[1,2] 
        𝔻𝒖ofw.value[index] = InverseJacobian[1,3] 

        𝔻𝒗ofu.value[index] = InverseJacobian[2,1] 
        𝔻𝒗ofv.value[index] = InverseJacobian[2,2] 
        𝔻𝒗ofw.value[index] = InverseJacobian[2,3] 

        𝔻𝐰ofu.value[index] = InverseJacobian[3,1] 
        𝔻𝐰ofv.value[index] = InverseJacobian[3,2] 
        𝔻𝐰ofw.value[index] = InverseJacobian[3,3] 

        DxofX = [𝔻uof𝒖.value[index] 𝔻uof𝒗.value[index] 𝔻uof𝐰.value[index];
                 𝔻vof𝒖.value[index] 𝔻vof𝒗.value[index] 𝔻vof𝐰.value[index];
                 𝔻wof𝒖.value[index] 𝔻wof𝒗.value[index] 𝔻wof𝐰.value[index]]

        DXofx = [𝔻𝒖ofu.value[index] 𝔻𝒖ofv.value[index] 𝔻𝒖ofw.value[index];
                 𝔻𝒗ofu.value[index] 𝔻𝒗ofv.value[index] 𝔻𝒗ofw.value[index]; 
                 𝔻𝐰ofu.value[index] 𝔻𝐰ofv.value[index] 𝔻𝐰ofw.value[index]]
    end

    𝔻𝒖    = 𝔻𝒖ofu * 𝔻u + 𝔻𝒖ofv * 𝔻v + 𝔻𝒖ofw * 𝔻w
    𝔻𝒗    = 𝔻𝒗ofu * 𝔻u + 𝔻𝒗ofv * 𝔻v + 𝔻𝒗ofw * 𝔻w
    𝔻𝐰    = 𝔻𝐰ofu * 𝔻u + 𝔻𝐰ofv * 𝔻v + 𝔻𝐰ofw * 𝔻w
    
    return(𝔻𝒗, 𝔻𝒖, 𝔻𝐰)
end
