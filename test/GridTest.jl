#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

function check_interpolatePatch(N::Int, M::Int, loc::Array{Int,1})::Float64
    fNpatch  = Float64[sin(pi*x) + sin(pi*y) for x in chebgrid(N, M, loc[1]), y in chebgrid(N, M, loc[2])]
    f2Npatch = Float64[sin(pi*x) + sin(pi*y) for x in chebgrid(2*N, M, loc[1]), y in chebgrid(2*N, M, loc[2])]
    patch    = Patch([1,1], fNpatch)
    interpfNpatch = interpolatePatch(patch, chebgrid(2*N), chebgrid(2*N))
    return L2norm(interpfNpatch.value, f2Npatch, chebweights(2*N))
end

@test_broken check_interpolatePatch(12, 1, [1,1]) < 1e-14
@test_broken check_interpolatePatch(12, 2, [1,1]) < 1e-14
@test_broken check_interpolatePatch(12, 2, [1,2]) < 1e-14
@test_broken check_interpolatePatch(12, 2, [2,1]) < 1e-14
@test_broken check_interpolatePatch(12, 2, [2,2]) < 1e-14
