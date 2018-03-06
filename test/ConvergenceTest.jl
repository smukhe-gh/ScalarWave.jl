#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 02-2018
#--------------------------------------------------------------------

# Test p-convergence on -1 to 1
# Choose N points > restrict poly onto grid > 
# compare L2 norm by evaluating at 2N points
function testpconv(fn::Function, Nx::Int, Ny::Int)::Float64
    fpatch  = Float64[fn(x,y) for x in chebgrid(2Nx), y in chebgrid(2Ny)]
    frpatch = Patch([1,1], projectonPatchbyRestriction(fn, Nx, Ny))
    return L2norm(fpatch, interpolatePatch(frpatch, 2Nx, 2Ny).value, chebweights(2Nx), chebweights(2Ny))
end

# Test p-convergence on -1 to 1
# Choose N points > restrict poly onto grid > 
# compare with exact solution at 2N points
function testpconv(fn::Function, Nx::Int, Ny::Int, M::Int, loc::Array{Int,1})::Float64
    fpatch  = Float64[fn(x,y) for x in chebgrid(2Nx, M, loc[1]), y in chebgrid(2Ny, M, loc[2])]
    frpatch = Patch([1,1], projectonPatchbyRestriction(fn, Nx, Ny, M, loc))
    return L2norm(fpatch, interpolatePatch(frpatch, 2Nx, 2Ny).value, chebweights(2Nx)/M, chebweights(2Ny)/M)
end

# Test h-convergence
# restrict poly onto each patch > compare with exact poly prolongated onto patches
# sum all L2 norms at 2N points in each patch

function testhconv(fn::Function, Nx::Int, Ny::Int, M::Int)::Float64 
    L2err  = 0.0
    for m in 1:M, n in 1:M
        rPatch = Patch([1,1], projectonPatchbyRestriction(fn, Nx, Ny, M, [m,n]))
        sPatch = Float64[fn(x,y) for x in chebgrid(2Nx, M, m), y in chebgrid(2Ny, M, n)]
        L2err += L2norm(sPatch, interpolatePatch(rPatch, 2Nx, 2Ny).value, chebweights(2Nx)/M, chebweights(2Ny)/M)
    end
    return L2err
end

print("------------------------------------------------\n")
print("p-convergence in 1 to -1\n")
print("------------------------------------------------\n")
for n in 1:12
    np = 1
    f(x,y) = x^9 + y^4 + x^3*y^7
    L2error = testpconv(f, n, n)
    @show n, np^2, L2error
end

#=
print("\n------------------------------------------------\n")
print("p-convergence on arbitrary patch\n")
print("------------------------------------------------\n")
olderror = 1
for np in 1:44
    n = 4
    f(x,y) = x^9 + y^4 + x^3*y^7
    L2error = testpconv(f, n, n, np, [1,1])
    @show n, np^2, L2error
end
=#

print("\n------------------------------------------------\n")
print("h-convergence\n")
print("------------------------------------------------\n")
for np in 1:12
    n = 4
    f(x,y) = x^9 + y^4 + x^3*y^7
    L2error = testhconv(f, n, n, np)
    @show n, np^2, L2error
end 
