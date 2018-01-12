#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

function extractBC{T<:Integer}(patch::Array{Float64, 1}, s::T, N::T)
	boundary = zeros((N+1)^2)
	for i in 1:(N+1), j in 1:(N+1)
		I = (N+1)*(i-1) + j
		if s==0 && i==1
			boundary[I] = patch[I]
		elseif s==1 && j==1
			boundary[I] = patch[I]
		else
			continue		
		end
	end
	return boundary
end

function setBC!{T<:Integer}(patch::Array{Float64, 1}, boundary::Array{Float64, 1}, s::T, N::T)
	for i in 1:(N+1), j in 1:(N+1)
		I = (N+1)*(i-1) + j
		if s==0 && i==1
			patch[I] = boundary[I]
		elseif s==1 && j==1
			patch[I] = boundary[I]
		else
			continue		
		end
    end
    return patch
end

