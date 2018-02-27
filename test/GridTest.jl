#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

function testprolongate(fn::Function, Nx::Int, Ny::Int, M::Int)::Float64
    ef1patch = Float64[fn(x,y) for x in chebgrid(Nx), y in chebgrid(Ny)]
    nfMpatch = prolongate(Patch([1,1], ef1patch), M)
    efMpatch = zeros((Nx+1)*M, (Ny+1)*M)
    L2       = Array{Float64,1}[]
    for m in 1:M, n in 1:M
        efMpatch[1+(m-1)*Nx:1+m*Nx, 1+(n-1)*Ny:1+n*Ny] = Float64[fn(x,y) for x in chebgrid(Nx,M,m), y in chebgrid(Ny,M,n)]
    end
    return LInfnorm(efMpatch, nfMpatch)
end

function testrestrict(fn::Function, Nx::Int, Ny::Int, M::Int)::Float64
    dbase = Dict{Array{Int,1}, Patch}()
    for m in 1:M, n in 1:M
        dbase[[m,n]] = Patch([m,n], Float64[fn(x,y) for x in chebgrid(Nx,M,m), y in chebgrid(Ny,M,n)])
    end
    ef1patch = Float64[fn(x,y) for x in chebgrid(Nx), y in chebgrid(Ny)]
    nf1patch = restrict(dbase, M)
    return LInfnorm(ef1patch, nf1patch)
end

@test_broken testprolongate((x,y)->x^2+y^3, 4, 6, 2) < 1e-14
@test_broken testrestrict((x,y)->x^2+y^3, 4, 6, 2) < 1e-14
