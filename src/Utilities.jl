#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

function coordtrans{T<:Integer}(point::Array{Float64, 1}, loc::Array{T, 1}, M::T)
	x, y = point
    s = Float64[d for d in M-1:-2:-M+1]
	xp = (x + s[loc[1]])/M
	yp = (x + s[loc[2]])/M
	return xp, yp
end

function computeB{T<:Integer}(N::T, M::T, loc::Array{T, 1})
	B = initializeB(N, loc, M)
	if sum(loc) == 2
		continue
	elseif loc[1] == 1 || loc[2] == 1
		if loc[1] > loc[2]
			setBC!(B, extractBC(dbase[loc-[1,0]]), 0), 0, N)	
		else
			setBC!(B, extractBC(dbase[loc-[0,1]]), 1), 1, N)
		end
	else
		setBC!(B, extractBC(dbase[loc-[1,0]]), 0), 0, N)
		setBC!(B, extractBC(dbase[loc-[0,1]]), 1), 1, N)
	end 
	return B
end