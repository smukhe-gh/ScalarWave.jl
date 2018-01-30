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
	        operator[index] = (2.0/M)*chebw(i,N)*chebw(k,N)*chebd(k,l,N)*chebd(i,j,N)	
        end
	end
	return operator
end

function getIC{T<:Integer}(N::T, M::T, loc::Array{T,1}, fn::Function, s::Symbol)::Boundary	
    if s==:R
        xp = Float64[coordtrans(M, [chebx(i,N),chebx(1,N)], loc)[1] for i in 1:N+1]
        return Boundary(:R, fn.(xp))
    elseif s==:C
        yp = Float64[coordtrans(M, [chebx(1,N),chebx(j,N)], loc)[2] for j in 1:N+1]
        return Boundary(:C, fn.(yp))
    else
        error("Unknown symbol passed.")
    end
end
