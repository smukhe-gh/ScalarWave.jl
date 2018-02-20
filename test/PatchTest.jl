#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

function check_calcPatch(N::Int, bnd1::Function, bnd2::Function)::Float64
    fpatch = Float64[bnd1(x) + bnd2(y) for x in chebgrid(N), y in chebgrid(N)]
    patch  = Patch([1,1], fpatch) 
    bndx   = Boundary(:R, Float64[bnd1(x) for x in chebgrid(N)])
    bndy   = Boundary(:C, Float64[bnd2(y) for y in chebgrid(N)])
    error  = calcPatch([1,1], bndx, bndy, operator(N,1)).value - fpatch
    @test patch.value == fpatch
    @test patch.loc   == [1,1]
    @test getPB(patch, :R).value == fpatch[end, :]
    @test getPB(patch, :C).value == fpatch[:, end]
    return  L2norm(calcPatch([1,1], bndx, bndy, operator(N,1)).value, fpatch, chebweights(N))
end

@test check_calcPatch(12, x->sin(pi*x), y->sin(pi*y)) < 1e-14
@test_broken check_calcPatch(12, x->x, y->y) < 1e-14

