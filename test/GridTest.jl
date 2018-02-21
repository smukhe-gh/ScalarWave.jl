#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

function check_interpolatePatch(N::Int, M::Int, loc::Array{Int,1})::Float64
    fNpatch  = Float64[x^2 + y^3 + x^3*y^2 for x in chebgrid(N, M, loc[1]), y in chebgrid(N, M, loc[2])]
    f2Npatch = Float64[x^2 + y^3 + x^3*y^2 for x in chebgrid(2N, M, loc[1]), y in chebgrid(2N, M, loc[2])]
    patch    = Patch([1,1], fNpatch)
    interpfNpatch = interpolatePatch(patch, 2N)
    return L2norm(interpfNpatch.value, f2Npatch, chebweights(2*N))
end

@test check_interpolatePatch(20, 1, [1,1]) < 1e-14
@test check_interpolatePatch(20, 2, [1,1]) < 1e-14
@test check_interpolatePatch(20, 2, [1,2]) < 1e-14
@test check_interpolatePatch(20, 2, [2,1]) < 1e-14
@test check_interpolatePatch(20, 2, [2,2]) < 1e-14
