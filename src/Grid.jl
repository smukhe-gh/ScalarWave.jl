#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
# XXX: WORK IN PROGRESS
#--------------------------------------------------------------------

include("Physics.jl")
include("Patch.jl")

# compute each patchindex in order and save the results in a dictionary
# XXX: Use spawn and fetch to implement futures
function distribute(N::Int64, M::Int64)

	grid  = Int64[i + j - 1 for i in 1:M, j in 1:M]
	A     = Physics.operator(N)	# call the operator beforehand.
	# Dict{CartesianIndex{Int, 2}, Any}
	dbase = Dict{Array{Int64, 1}, Any}()	# Set type of remote references.

	for i in 1:2M-1
		slice = map(x->collect(ind2sub(grid, x)), find(grid .== i))
	
		for patchindex in slice
			b = zeros(N+1, N+1)	# boundary vector (passed by reference; get's changed)
			initialize(b)		# this sets all the background potential and BCs.
								# b needs to be a vector from start.

			# Also, have a higher layer, which treats everything in N-D and
			# re-interpret it has lower level strucutres (2D for the operator)
			# and vector for b.
			if sum(patchindex) == 2
				continue
			elseif patchindex[1] == 1 || patchindex[2] == 1
				if patchindex[1] > patchindex[2]
					# FIXME: Not the correct implementation! 
					# Currently it will stop here and wait for the call to finish.
					Patch.setBC(b, Patch.extractBC(fetch(dbase[patchindex-[1,0]]), 0), 0)	
				else
					Patch.setBC(b, Patch.extractBC(fetch(dbase[patchindex-[0,1]]), 1), 1)
				end
			else
				Patch.setBC(b, Patch.extractBC(fetch(dbase[patchindex-[1,0]]), 0), 0)
				Patch.setBC(b, Patch.extractBC(fetch(dbase[patchindex-[0,1]]), 1), 1)
			end 

			# store the future of the result in the dbase. 
			# XXX: vec(b') preserves the correct order of raveling.
			dbase[patchindex] = @spawn A\vec(b')
		end
	end
	return tree
end


# Split this routine into several concerns. 
# 	- Extract (set b in one function)
# 	- Solve 
#   - Parallelization (serial and parallel wrappers/functors)


