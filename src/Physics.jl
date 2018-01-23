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

function getIC{T<:Integer}(N::T, M::T, loc::Array{T,1}, s::Symbol)::Boundary	
	B = zeros(N+1)
    fnbrow(x,y) = sin(pi*x)
    fnbcol(x,y) = sin(pi*y)
    for i in 1:N+1
        if s==:R
            xp, yp = coordtrans(M, [chebx(i,N),chebx(1,N)], loc)
			B[i]   = fnbrow(xp,yp) 
        else
            xp, yp = coordtrans(M, [chebx(1,N),chebx(i,N)], loc)
			B[i]   = fnbcol(xp,yp) 
		end   
	end
    if s==:R
	    return Boundary(:R, B)
    else
        return Boundary(:C, B)
    end
end
