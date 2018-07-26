#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 02-2018
#--------------------------------------------------------------------

function testdrawmultipatch()
    P = 20
    dict = distribute(u->0, v->exp(-v^2/0.1),
                     (u,v)-> 0,
                     P, P, 1)
    
    patch = dict[[1,1]]
    u = chebgrid(P)
    plot(u, patch.value[1,:], "m-")
    plot(u, patch.value[end,:],"ro")
    savefig("boundary-profiles.pdf")
    close()
    drawmultipatch(dict, "visualization-test")
    return true
end

@test testdrawmultipatch() == true
#run(`rm visualization-test.pdf`)

