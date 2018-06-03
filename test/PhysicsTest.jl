#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

function testoperator(Nx::Int, Ny::Int)::Array{Float64,4}
    """
    The main operator we're using 
    sum(chebw(i,Nx)*chebw(j,Ny)*
         chebd(i,m,Nx)*delta(j,n)*phi[m,n]*
         delta(i,p)*chebd(j,q,Ny)*phi[p,q] 
         for m in 1:Nx+1, n in 1:Ny+1, p in 1:Nx+1, q in 1:Ny+1, 
         i in 1:Nx+1, j in 1:Ny+1)
    where we sum over i, j in the computation of the operator.
    
    And for computing the action of the operator on the vector field
    operator[m,n,p,q] phi[m,n] phi[p,q]
    and when you reshape it into a 2D array (m,n) and (p,q) are bunched 
    together.
    """

    operator = zeros(Nx+1, Ny+1, Nx+1, Ny+1)
    for index in CartesianRange(size(operator)) 
        (m,n,p,q) = index.I
        operator[m,n,p,q] = chebw(p,Nx)*chebw(n,Ny)*chebd(p,m,Nx)*chebd(n,q,Ny)
    end
    return operator
end

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

@test_broken 2*testoperator(2, 2) ≈ derivOP(2, 2)
@test_broken 2*testoperator(2, 4) ≈ derivOP(2, 4)
@test testboundaryOP(2,2) == shapeH2L(boundaryOP(2,2))
@test testboundaryOP(2,4) == shapeH2L(boundaryOP(2,4))
@test RHS((x,y)->x^8+y^9, 8, 4, 7, [2,1]) ≈ testRHS((x,y)->x^8+y^9, 8, 4, 7, [2,1])
