#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

function operator{T<:Int}(N::T, M::T)::Array{Float64, 4}
	operator = zeros(N+1, N+1, N+1, N+1)
    for index in CartesianRange(size(operator))    
        k = index.I[1]
        i = index.I[2]
        l = index.I[3]
        j = index.I[4]   
        if 	i==1 || k==1
		    operator[index] = delta(i,j)*delta(k,l)
        else 
	        operator[index] = (2.0/M^4)*chebw(i,N)*chebw(k,N)*chebd(k,l,N)*chebd(i,j,N)	
        end
	end
	return operator
end

function getIC{T<:Integer}(N::T, M::T, loc::Array{T,1}, fn::Function, s::Symbol)::Boundary	
    if s==:R
        xg = chebgrid(N, M, loc[1]) 
        return Boundary(:R, fn.(xg))
    elseif s==:C
        yg = chebgrid(N, M, loc[2]) 
        return Boundary(:C, fn.(yg))
    else
        error("Unknown symbol passed.")
    end
end

