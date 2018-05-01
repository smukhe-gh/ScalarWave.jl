#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 02-2018
#--------------------------------------------------------------------

function testdrawmultipatch()
    dict = distribute(u->0, v->0, (u,v)->exp(-u^2/0.1)*exp(-v^2/0.1), 4, 4, 2)
    drawmultipatch(dict)
    return true
end

@test testdrawmultipatch() == true
