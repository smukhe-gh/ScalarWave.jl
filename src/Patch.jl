#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

function extractBC{T<:Integer}(patch::Array{Float64,2}, s::T)::Array{Float64,1}
    N = size(patch)[1] - 1
    boundary = zeros(N+1)
	for i in 1:(N+1), j in 1:(N+1)
		if s==1 && i==1
			boundary[j] = patch[i,j]
		elseif s==0 && j==1
			boundary[i] = patch[i,j]
		else
			continue		
		end
	end
	return boundary
end

function setBC!{T<:Integer}(patch::Array{Float64,2}, boundary::Array{Float64,1}, s::T)::Array{Float64,2}
	N = size(patch)[1] - 1
    for i in 1:(N+1), j in 1:(N+1)
		if s==0 && i==1
			patch[i,j] = boundary[j]
		elseif s==1 && j==1
			patch[i,j] = boundary[i]
		else
			continue		
		end
    end
    return patch
end

