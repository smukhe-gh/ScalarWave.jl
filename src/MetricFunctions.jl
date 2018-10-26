#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2018
# Define operations on the metric 
#--------------------------------------------------------------------

import LinearAlgebra: inv, det, eigvals

function inv(g::Metric{_dd, 4})::Metric{_uu, 4}   
    ginv = Metric{_uu, 4}([similar(g[1,1]),
                          similar(g[2,1]), similar(g[2,2]),
                          similar(g[3,1]), similar(g[3,2]), similar(g[3,3]),
                          similar(g[4,1]), similar(g[4,2]), similar(g[4,3]), similar(g[4,4])])

    for index in CartesianIndices(size(g[1,1].space))
        # TODO: Remove array allocations from inside the loop to make things faster.
        tempg = [g[1,1].value[index]  g[1,2].value[index] g[1,3].value[index] g[1,4].value[index];
                 g[2,1].value[index]  g[2,2].value[index] g[2,3].value[index] g[2,4].value[index];
                 g[3,1].value[index]  g[3,2].value[index] g[3,3].value[index] g[3,4].value[index];
                 g[4,1].value[index]  g[4,2].value[index] g[4,3].value[index] g[4,4].value[index]]

        tempginv = inv(tempg)    

        ginv[1,1].value[index] = tempginv[1,1] 
        ginv[2,1].value[index] = tempginv[2,1] 
        ginv[3,1].value[index] = tempginv[3,1] 
        ginv[4,1].value[index] = tempginv[4,1] 
                                             
        ginv[2,2].value[index] = tempginv[2,2] 
        ginv[3,2].value[index] = tempginv[3,2] 
        ginv[4,2].value[index] = tempginv[4,2] 
                                             
        ginv[3,3].value[index] = tempginv[3,3] 
        ginv[4,3].value[index] = tempginv[4,3] 
                                             
        ginv[4,4].value[index] = tempginv[4,4] 

    end
    return ginv
end

function det(g::Metric{_dd, 4})::Field
    detg = similar(g[1,1])

    for index in CartesianIndices(size(g[1,1].space))
        # TODO: Remove array allocations from inside the loop to make things faster.
        tempg = [g[1,1].value[index]  g[1,2].value[index] g[1,3].value[index] g[1,4].value[index];
                 g[2,1].value[index]  g[2,2].value[index] g[2,3].value[index] g[2,4].value[index];
                 g[3,1].value[index]  g[3,2].value[index] g[3,3].value[index] g[3,4].value[index];
                 g[4,1].value[index]  g[4,2].value[index] g[4,3].value[index] g[4,4].value[index]]
        detg.value[index] = det(tempg)    
    end
    return detg
end

function eigvals(g::Metric{_dd, 4})
    E1 = similar(g[1,1])
    E2 = similar(g[1,1])
    E3 = similar(g[1,1])
    E4 = similar(g[1,1])
    
    for index in CartesianIndices(size(g[1,1].space))
        # TODO: Remove array allocations from inside the loop to make things faster.
        tempg = [g[1,1].value[index]  g[1,2].value[index] g[1,3].value[index] g[1,4].value[index];
                 g[2,1].value[index]  g[2,2].value[index] g[2,3].value[index] g[2,4].value[index];
                 g[3,1].value[index]  g[3,2].value[index] g[3,3].value[index] g[3,4].value[index];
                 g[4,1].value[index]  g[4,2].value[index] g[4,3].value[index] g[4,4].value[index]]
        E1, E2, E3, E4 = eigvals(tempg)    
    end

    return (E1, E2, E3, E4)
end
