#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 02-2018
#--------------------------------------------------------------------

function testpconv(fn::Function, Nx::Int, Ny::Int)::Float64
    fpatch  = Float64[fn(x,y) for x in chebgrid(2Nx), y in chebgrid(2Ny)]
    frpatch = Patch([1,1], projectonPatchbyRestriction(fn, Nx, Ny))
    return L2norm(fpatch, interpolatePatch(frpatch, 2Nx, 2Ny).value, chebweights(2Nx), chebweights(2Ny))
end

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
for n in 1:24
    np = 1
    f(x,y) = x^9 + y^4 + x^3*y^7
    L2error = testpconv(f, n, n)
    @show n, np^2, L2error
end

print("\n------------------------------------------------\n")
print("h-convergence\n")
print("------------------------------------------------\n")

L20 = 1
for np in 1:5
    n = 2
    f(x,y) = x^9 + y^4 + x^3*y^7
    L2n = testhconv(f, n, n, 3^np)
    @show n, (3^np)^2, L2n, L20/L2n
    L20 = L2n
end 
