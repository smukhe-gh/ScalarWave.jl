#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

function distribute{T<:Integer}(N::T, M::T)
	domain  = Int64[i + j - 1 for i in 1:M, j in 1:M]
	A     = operator(N)					# compute operator
	dbase = Dict{Array{Int, 1}, Array{Float64, 1}}()
	
	for i in 1:2M-1
		slice = map(x->collect(ind2sub(domain, x)), find(domain .== i))	
		for loc in slice
			B = computeB(N, loc, M)		# compute b (Ax=b)
			dbase[loc] = gesv!(A, B)	# solve the system (Note: A is overwritten with 
										# its LU factorization and B is overwritten with the solution X)
		end
	end	
	return dbase
end




