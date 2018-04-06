#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

function fdistribute{T<:Integer}(fbndr::Function, fbndc::Function, frhs::Function, Nx::T, Ny::T, M::T)::Dict{Array{Int,1}, Patch}
    dop = derivOP(Nx, Ny)
    bop = boundaryOP(Nx, Ny)
    fdbase = Dict{Array{Int, 1}, Future}()
    dbase  = Dict{Array{Int, 1}, Patch}()

    # This should run in parallel, scattered over multiple processes.
    for i in 2:2M, k in i-min(i-1,M):min(i-1,M) 
        loc  = [k, i-k]
        rhs  = fRHS(frhs, Nx, Ny, M, loc)
        bndx = (loc[2]==1) ? fgetPatchIC(fbndr, 0, Nx, M, loc[1]) : fgetPatchBnd(fdbase[loc-[0,1]], 0)
        bndy = (loc[1]==1) ? fgetPatchIC(fbndc, 1, Ny, M, loc[2]) : fgetPatchBnd(fdbase[loc-[1,0]], 1)
        fdbase[loc] = fcalcPatch(bndx, bndy, rhs, dop, bop, loc)
    end
    
    # This is where we collect all the results
    # We could improve this. Instead of fetching all, we could fetch only the last one.
    # since the previous ones are already fetched.
    for i in 2:2M, k in i-min(i-1,M):min(i-1,M) 
        loc  = [k, i-k]
        dbase[loc] = fetch(fdbase[loc]) 
    end

    return dbase
end

function fgetPatchIC{T<:Integer}(fbndr::Function, s::T, N::T, M::T, loc::T)::Future
    @spawn getPatchIC(fbndr, s, N, M, loc)
end

function fRHS{T<:Integer}(frhs::Function, Nx::T, Ny::T, M::T, loc::Array{Int,1})::Future
    @spawn RHS(frhs, Nx, Ny, M, loc)
end

function fgetPatchBnd(fpatch::Future, s::Int)::Future
    @spawn getPatchBnd(fetch(fpatch), s)
end

function fcalcPatch{T<:Array{Float64,4}}(fbndx::Future, fbndy::Future, frhs::Future, dop::T, bop::T, loc::Array{Int,1})::Future
    @spawn calcPatch(fetch(fbndx), fetch(fbndy), fetch(frhs), dop, bop, loc)
end
