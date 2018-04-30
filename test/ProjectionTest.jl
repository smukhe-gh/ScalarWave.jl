#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

function testprojectonPatchBndbyRestriction(fn::Function, Nx::Int)::Float64
    pfbnd = projectonPatchBndbyRestriction(fn, Nx) 
    efbnd = Float64[fn(x) for x in chebgrid(Nx)]
    w     = chebweights(Nx)  
    L2error = sqrt(w'*(pfbnd-efbnd).^2)
    return L2error 
end

function testprojectonPatchBndbyRestriction(fn::Function, Nx::Int, M::Int, loc::Int)::Float64
    pfbnd = projectonPatchBndbyRestriction(fn, Nx, M, loc) 
    efbnd = Float64[fn(x) for x in chebgrid(Nx, M, loc)]
    w     = chebweights(Nx)  
    L2error = sqrt(w'*(pfbnd-efbnd).^2)
    return L2error 
end

function testprojectonPatchbyRestriction(fn::Function, Nx::Int, Ny::Int)::Float64
    pfpatch = projectonPatchbyRestriction(fn, Nx, Ny)
    efpatch = Float64[fn(x,y) for x in chebgrid(Nx), y in chebgrid(Ny)]
    return L2norm(pfpatch, efpatch, chebweights(Nx), chebweights(Ny))
end

function testprojectonPatchbyRestriction(fn::Function, Nx::Int, Ny::Int, M::Int, loc::Array{Int,1})::Float64
    pfpatch = projectonPatchbyRestriction(fn, Nx, Ny, M, loc)
    efpatch = Float64[fn(x,y) for x in chebgrid(Nx, M, loc[1]), y in chebgrid(Ny, M, loc[2])]
    return L2norm(pfpatch, efpatch, chebweights(Nx), chebweights(Ny))
end

@test testprojectonPatchBndbyRestriction(x->x^8, 7) > testprojectonPatchBndbyRestriction(x->x^8, 8)
@test testprojectonPatchBndbyRestriction(x->x^8, 6, 4, 1) >  testprojectonPatchBndbyRestriction(x->x^8, 6, 8, 1)
@test testprojectonPatchBndbyRestriction(x->x^8, 7, 4, 1) >  testprojectonPatchBndbyRestriction(x->x^8, 8, 4, 1)
@test testprojectonPatchbyRestriction((x,y)->x^8+y^9, 8, 8) > testprojectonPatchbyRestriction((x,y)->x^8+y^9, 8, 9)
@test testprojectonPatchbyRestriction((x,y)->x^8+y^9, 8, 4, 7, [2,1]) > testprojectonPatchbyRestriction((x,y)->x^8+y^9, 5, 9, 3, [2,1])
@test projectonPatchbyRestriction((x,y)->x+y, 2, 2, 2, [1,1])[:, end] ≈ projectonPatchbyRestriction((x,y)->x+y, 2, 2, 2, [1,2])[:, 1]
@test projectonPatchbyRestriction((x,y)->x+y, 2, 2, 2, [1,1])[end, :] ≈ projectonPatchbyRestriction((x,y)->x+y, 2, 2, 2, [2,1])[1, :]

