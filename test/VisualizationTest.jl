#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 02-2018
#--------------------------------------------------------------------

function testdrawpatch()
    drawpatch(Float64[x+y for x in chebgrid(20), y in chebgrid(20)])
    run(`rm ./luxor-patch.pdf`)
    return true
end

@test testdrawpatch() == true
