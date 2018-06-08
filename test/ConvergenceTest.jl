#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 02-2018
#--------------------------------------------------------------------
"TODO: Fix convergence for n > 243
       Test convergence for analytic potential
       Test self-convergence with arbitrary potentials"

function testpconvforprojection(fn::Function, Nx::Int, Ny::Int)::Float64
    rPatch = Patch([1,1], projectonPatchbyRestriction(fn, Nx, Ny, 1, [1,1]))
    sPatch  = Float64[fn(x,y) for x in chebgrid(4Nx), y in chebgrid(4Ny)]
    return L2norm(sPatch, interpolatePatch(rPatch, 4Nx, 4Ny).value, chebweights(4Nx), chebweights(4Ny))
end

function testhconvforprojection(fn::Function, Nx::Int, Ny::Int, M::Int)::Float64 
    L2err  = 0.0
    for m in 1:M, n in 1:M
        rPatch = Patch([1,1], projectonPatchbyRestriction(fn, Nx, Ny, M, [m,n]))
        sPatch = Float64[fn(x,y) for x in chebgrid(4Nx, M, m), y in chebgrid(4Ny, M, n)]
        L2err += (L2norm(sPatch, interpolatePatch(rPatch, 4Nx, 4Ny).value, chebweights(4Nx)/M, chebweights(4Ny)/M))^2
    end
    return sqrt(L2err)
end

# test p-refinement
showconv(x->sin(pi*x), y->sin(pi*y), 
        (x,y)-> 0, 
        (x,y)-> sin(pi*x) + sin(pi*y),
        20,  # maxmodes 
        1,   # maxlevels                      
        1)   # h-factor                               

# test h-refinement
showconv(x->sin(pi*x), y->sin(pi*y), 
        (x,y)-> 0, 
        (x,y)-> sin(pi*x) + sin(pi*y),
        2,  # maxmodes 
        8,  # maxlevels                      
        2)  # h-factor                               
        
