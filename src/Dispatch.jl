#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

function distribute{T<:Integer}(fbndr::Function, fbndc::Function, frhs::Function, Nx::T, Ny::T, M::T)::Dict
	dop = derivOP(Nx, Ny)/M^4
    bop = boundaryOP(Nx, Ny)
    dbase  = Dict{Array{Int, 1}, Patch}()
    for i in 2:2M, k in i-min(i-1,M):min(i-1,M) 
        loc  = [k, i-k]
        rhs  = RHS(frhs, Nx, Ny, M, loc)
        brow = (loc[1]==1) ? getIC(fbndr, Nx, M, loc[1]) : getPB(dbase[loc-[1,0]], 0)
        bcol = (loc[2]==1) ? getIC(fbndc, Ny, M, loc[2]) : getPB(dbase[loc-[0,1]], 1)
        dbase[loc] = calcPatch(loc, brow, bcol, op)
	end
    return dbase
end
