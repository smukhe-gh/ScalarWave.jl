#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

type BC = Array{Float64, 1}
type RHS = Array{Float64, 1}
# function computeB(N, loc, M, BC):RHS
# function extractBC(SOL):BC

function fextractBC(fsol::RemoteRef{SOL}): RemoteRef{BC}
	@spawn extractBC(fetch(fsol))
end

function fcomputeB(N, loc, M, fbc::RemoteRef{BC}): RemoteRef{RHS}
	@spawn computeB(N, loc, M, fetch(fbc))
end

function distribute{T<:Integer}(N::T, M::T)
	domain  = Int64[i + j - 1 for i in 1:M, j in 1:M]
	A     = operator(N)					# compute operator
	# Dict: Array{Int} -> CartesianIndex
	dbase = Dict{Array{Int, 1}, Array{Float64, 1}}()
	
	for i in 1:2M-1
		# use filter instead
		slice = map(x->collect(ind2sub(domain, x)), find(domain .== i))	
		for loc in slice
			B = computeB(N, loc, M)		# compute b (Ax=b)
			dbase[loc] = gesv!(A, B)	# solve the system (Note: A is overwritten with 
										# its LU factorization and B is overwritten with the solution X)
		end
	end	
	return dbase
end

function distribute2(M)
	allinds = Tuple{Int,Int}[]
	for sij in 1:2M+1
		# do stuff
		for dij in -sij:sij
			(i2, j2) = (sij+dij, sij-dij)
			if i2 % 2 == 0 && j2 % 2 == 0
				(i, j) = (i `div` 2, j `div` 2)
				if i>=1 && i<=M && j>=1 && j<=M
					# do per-patch stuff
					push!(allinds, (i,j))
				end
			end
		end
	end
	sort!(allinds)
	println allinds
end


