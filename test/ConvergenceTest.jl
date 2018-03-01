#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 02-2018
#--------------------------------------------------------------------

# Test p-convergence on -1 to 1
# Choose N points > restrict poly onto grid > 
# compare with exact solution at N points
function testpconv(fn::Function, Nx::Int, Ny::Int)::Float64
    fpatch  = Float64[fn(x,y) for x in chebgrid(Nx), y in chebgrid(Ny)]
    frpatch = projectonPatchbyRestriction(fn, Nx, Ny)
    return L2norm(fpatch, frpatch, chebweights(Nx), chebweights(Ny))
end

print("------------------------------------------------\n")
print("p-convergence in 1 to -1\n")
print("------------------------------------------------\n")
for n in 2:12
    f(x,y) = x^9 + y^4 + x^3*y^7
    L2error = testpconv(f, n, n)
    @show n, L2error
end

# Test p-convergence on -1 to 1
# Choose N points > restrict poly onto grid > 
# compare with exact solution at N points
function testpconv(fn::Function, Nx::Int, Ny::Int, M::Int, loc::Array{Int,1})::Float64
    fpatch  = Float64[fn(x,y) for x in chebgrid(Nx, M, loc[1]), y in chebgrid(Ny, M, loc[2])]
    frpatch = projectonPatchbyRestriction(fn, Nx, Ny, M, loc)
    return L2norm(fpatch, frpatch, chebweights(Nx)/M, chebweights(Ny)/M)
end

print("\n------------------------------------------------\n")
print("p-convergence on arbitrary patch\n")
print("------------------------------------------------\n")
olderror = 1
for m in 2:12
    f(x,y) = sin(x)^9 + y^4 + x^3*y^7
    L2error = testpconv(f, 4, 4, m, [2,1])
    @show m, L2error #, olderror/L2error
    olderror = L2error 
end


# Test h-convergence
# Choose N points > compute exact solution at N points > prolongate into many patches
# restrict poly onto each patch > compare with exact poly prolongated onto patches
# sum all L2 norms

function testhconv(fn::Function, Nx::Int, Ny::Int, M::Int)::Float64 
    sPatch = Float64[fn(x,y) for x in chebgrid(Nx), y in chebgrid(Ny)]
    mPatch = prolongatePatch(Patch([1,1], sPatch), M)
    L2err  = 0.0
    for m in 1:M, n in 1:M
        rPatch = projectonPatchbyRestriction(fn, Nx, Ny, M, [m,n])
        L2err += L2norm(mPatch[[m,n]].value, rPatch, chebweights(Nx)/M, chebweights(Ny)/M)
    end
    return L2err
end

print("\n------------------------------------------------\n")
print("h-convergence\n")
print("------------------------------------------------\n")
for m in 2:12
    f(x,y) = x^9 + y^4 + x^3*y^7
    L2error = testhconv(f, 7, 7, 2*m)
    @show 2*m, L2error
end 


