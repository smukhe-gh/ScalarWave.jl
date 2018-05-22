#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 02-2018
#--------------------------------------------------------------------

function testdrawmultipatch()
    dict = distribute(u->0, v->0, 
                     (u,v)-> -exp(-u^2 - v^2)*(4v*(u*cos(2u)-u*cos(2v)+sin(2u))- 4u*sin(2v)),
                     40, 40, 12)
    drawmultipatch(dict)
    return true
end

@test testdrawmultipatch() == true
