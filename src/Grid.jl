#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------
# TODO: Use Dict{Array{Int,1}, CartesianIndex}

function distribute{T<:Integer}(N::T, M::T)::Dict
	A = reshapeA(operator(N))
	dbase = Dict{Array{Int, 1}, Array{Float64, 1}}()
    for i in 2:2M, k in i-min(i-1,M):min(i-1,M) 
        loc = [k, i-k]	
		B   = computeRHS(N, M, loc, dbase)
		dbase[loc] = shapeB(solve(A, reshapeB(B)))
	end	
	return dbase
end
