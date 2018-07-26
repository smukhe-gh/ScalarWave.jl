#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

function testboundaryOP(Nx::Int, Ny::Int)::Array{Float64,2}
    bnd = zeros(Nx+1, Ny+1)
    bnd[1,:] = bnd[:,1] = 1
    return diagm(vec(bnd))
end    

function testRHS(fn::Function, Nx::Int, Ny::Int, M::Int, loc::Array{Int,1})
    patch = projectonPatchbyRestriction(fn, Nx, Ny, M, loc)
    for index in CartesianRange(size(patch))
        i = index.I[1]
        j = index.I[2]
        patch[index] = chebw(i,Nx)*chebw(j,Ny)*patch[index]
    end
    return patch
end

@test testboundaryOP(2,2) == shapeH2L(boundaryOP(2,2))
@test testboundaryOP(2,4) == shapeH2L(boundaryOP(2,4))
@test RHS((x,y)->x^8+y^9, 8, 4, 7, [2,1]) ≈ testRHS((x,y)->x^8+y^9, 8, 4, 7, [2,1])
