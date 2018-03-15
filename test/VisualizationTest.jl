#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 02-2018
#--------------------------------------------------------------------

function testdrawpatch()
    drawpatch(Float64[x+y for x in chebgrid(20), y in chebgrid(20)])
    run(`rm ./luxor-patch.pdf`)
    return true
end

function testdrawarray()
    drawarray(Float64[x+y for x in collect(linspace(0,100,100)), y in collect(linspace(0,100,100))])
    run(`rm ./luxor-array.pdf`)
    return true
end

function testdrawgrid()
    fn(x,y) = sin(pi*x) + sin(pi*y)
    dbase = Dict{Array{Int,1}, Patch}()
    for m in 1:2, n in 1:2
        dbase[[m,n]] = Patch([m,n], projectonPatchbyRestriction(fn, 12, 12, 2, [m, n]))
    end
    drawgrid(dbase)
    return true
end

@test testdrawpatch() == true
@test testdrawarray() == true
@test testdrawgrid()  == true
