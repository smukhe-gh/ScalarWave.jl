#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

function distribute{T<:Integer}(fbndr::Function, fbndc::Function, frhs::Function, Nx::T, Ny::T, M::T)::Dict{Array{Int,1}, Patch}
	dop = derivOP(Nx, Ny)/M^4
    bop = boundaryOP(Nx, Ny)
    dbase  = Dict{Array{Int, 1}, Patch}()
    for i in 2:2M, k in i-min(i-1,M):min(i-1,M) 
        loc  = [k, i-k]
        rhs  = RHS(frhs, Nx, Ny, M, loc)
        bndx = (loc[2]==1) ? getPatchIC(fbndr, 0, Nx, M, loc[1]) : getPatchBnd(dbase[loc-[0,1]], 0)
        bndy = (loc[1]==1) ? getPatchIC(fbndc, 1, Ny, M, loc[2]) : getPatchBnd(dbase[loc-[1,0]], 1)
        dbase[loc] = calcPatch(bndx, bndy, rhs, dop, bop, loc)
    end
    return dbase
end
