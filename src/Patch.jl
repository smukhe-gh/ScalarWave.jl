#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

module Patch
export extractBC, setBC

	function extractBC(patch::Array{Float64, 2}, s::Int64)
		s==0 ? return Float64[x for x in patch[end,:]] : return Float64[x for x in patch[:, end]]
	end

	function setBC(patch::Array{Float64, 2}, boundary::Array{Float64, 2}, s::Int64)
		if s==0 # set row
			for (index, value) in enumerate(boundary):
				patch[1, index] = value
		else	# set column
			for (index, value) in enumerate(boundary):
				patch[index, 1] = value
	end

	function computepatchcoefficents(patch::Array{Float64, 2})
		# see Mathematica notebook for index implementation
	end

end 	# end module