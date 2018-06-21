#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 02-2018
#--------------------------------------------------------------------

function testdrawmultipatch()
    dict = distribute(u->0*exp(-u^2/0.1), v->exp(-v^2/0.1),
                     (u,v)-> 0,
                     30, 30, 1)
    
    patch = dict[[1,1]]
    u = chebgrid(30)
    plot(u, patch.value[1,:])
    plot(u, patch.value[end,:])
    savefig("boundary-profiles.pdf")
    close()
    drawmultipatch(dict, "visualization-test")
    return true
end

@test testdrawmultipatch() == true
#run(`rm visualization-test.pdf`)

