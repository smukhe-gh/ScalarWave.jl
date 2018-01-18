#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

function extractBC(patch::Array{Float64,2}, s::Symbol)::Array{Float64,1}
    N = size(patch)[1] - 1
    boundary = zeros(N+1)
	for i in 1:(N+1), j in 1:(N+1)
		if s==:R && i==N+1
			boundary[j] = patch[i,j]
		elseif s==:C && j==N+1
			boundary[i] = patch[i,j]
		else
			continue		
		end
	end
	return boundary
end

function setBC!(patch::Array{Float64,2}, boundary::Array{Float64,1}, s::Symbol)::Array{Float64,2}
	N = size(patch)[1] - 1
    for i in 1:(N+1), j in 1:(N+1)
		if s==:R && i==1
			patch[i,j] = boundary[j]
		elseif s==:C && j==1
			patch[i,j] = boundary[i]
		else
			continue		
		end
    end
    return patch
end

function solve(A::Array{Float64,2}, RHS::Array{Float64,1})::Array{Float64,1}
    return A/RHS
end

