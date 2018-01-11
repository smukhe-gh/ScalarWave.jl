#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

module Utilities
export imat, affine

	function imat(i::Int64, j::Int64)
	    i==j ? 1.0 : 0.0
	end

	function affine(x::Float64, pindex::Int64, M::Int64)
		s = Float64[d for d in M-1:-2:-M+1]
		x = (x + s[pindex[1]])/M
		return x
	end

end 	# end module
