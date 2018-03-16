#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 02-2018
#--------------------------------------------------------------------

function testpconv(fn::Function, Nx::Int, Ny::Int)::Float64
    rPatch = Patch([1,1], projectonPatchbyRestriction(fn, Nx, Ny, 1, [1,1]))
    sPatch  = Float64[fn(x,y) for x in chebgrid(4Nx), y in chebgrid(4Ny)]
    return L2norm(sPatch, interpolatePatch(rPatch, 4Nx, 4Ny).value, chebweights(4Nx), chebweights(4Ny))
end

function testhconv(fn::Function, Nx::Int, Ny::Int, M::Int)::Float64 
    L2err  = 0.0
    for m in 1:M, n in 1:M
        rPatch = Patch([1,1], projectonPatchbyRestriction(fn, Nx, Ny, M, [m,n]))
        sPatch = Float64[fn(x,y) for x in chebgrid(4Nx, M, m), y in chebgrid(4Ny, M, n)]
        L2err += (L2norm(sPatch, interpolatePatch(rPatch, 4Nx, 4Ny).value, chebweights(4Nx)/M, chebweights(4Ny)/M))^2
    end
    return sqrt(L2err)
end

function testdistribute(bnd1::Function, bnd2::Function, Nx::Int, Ny::Int, M::Int)::Float64
    dbase = distribute(bnd1, bnd2, (x,y)->0, Nx, Ny, M)
    sumL2 = 0.0
    for m in 1:M, n in 1:M
        sPatch = interpolatePatch(dbase[[m,n]], 2Nx, 2Ny).value
        fPatch = Float64[bnd1(x) + bnd2(y) for x in chebgrid(2Nx, M, m), y in chebgrid(2Ny, M, n)]
        sumL2 += (L2norm(fPatch, sPatch, chebweights(2Nx)/M, chebweights(2Ny)/M))^2
    end
    return sqrt(sumL2)
end

#=
println("f(x,y) = x^9 + y^9")
print("------------------------------------------------\n")
print("p-convergence in 1 to -1\n")
print("------------------------------------------------\n")
for p in 1:10
    np = 1
    f(x,y) = x^9 + y^9
    L2error = testpconv(f, p, p)
    @show p, np, L2error
end

print("\n------------------------------------------------\n")
print("h-convergence\n")
print("------------------------------------------------\n")

L20 = 1
for np in 0:6
    p = 5
    f(x,y) = x^9 + y^9
    L2n = testhconv(f, p, p, 2^np)
    n = p + 1
    @show p, 2^np, L2n, L20/L2n
    L20 = L2n
end

print("\n------------------------------------------------\n")
print("p-convergence\n")
print("------------------------------------------------\n")

println("x->sin(pi*x), y->sin(pi*y)")
L20 = 1
for p in 0:12
    np  = 0
    p   = 2p + 1
    L2n = testdistribute(x->sin(pi*x), y->sin(pi*y), p, p, 3^np)
    @show p, 3^np, L2n
    L20 = L2n
end

print("\n------------------------------------------------\n")
print("h-convergence\n")
print("------------------------------------------------\n")

println("x->sin(pi*x), y->sin(pi*y)")
L20 = 1
for np in 0:4
    p = 3
    L2n = testdistribute(x->sin(pi*x), y->sin(pi*y), p, p, 3^np)
    @show p, 3^np, L2n, L20/L2n
    L20 = L2n
end
=#
