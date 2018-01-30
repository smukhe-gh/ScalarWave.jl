#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

function distribute{T<:Integer}(N::T, M::T, fr::Function, fc::Function)::Dict
	op = operator(N, M)
	dbase = Dict{Array{Int, 1}, Patch}()
    for i in 2:2M, k in i-min(i-1,M):min(i-1,M) 
        loc  = [k, i-k]	
        brow = (loc[1]==1) ? getIC(N, M, loc, fr, :R) : getPB(dbase[loc-[1,0]], :R)
        bcol = (loc[2]==1) ? getIC(N, M, loc, fc, :C) : getPB(dbase[loc-[0,1]], :C)
        dbase[loc] = calcPatch(loc, brow, bcol, op)
	end	
	return dbase
end
