#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 02-2018
#--------------------------------------------------------------------

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

#@test sconv(x->x^3-1, y->y^3-1, 3, 3, 2) < 1e-14
#showconv(x->x^5-1, y->y^5-1, 12, 1, 1)
#showconv(x->sin(pi*x), y->sin(pi*y), 2, 8, 2)
#showconv(x->sin(pi*x), y->sin(pi*y), 2, 6, 3)
showconv(v->0, v->0, 
        (u,v)-> exp(-u^2 - v^2)*(4v*(u*cos(2u)-u*cos(2v)+sin(2u))- 4u*sin(2v)),
        (v,u)-> sin(v-u)*sin(v+u)*exp(-(u^2 + v^2)), 
        3, 4, 3)
        
