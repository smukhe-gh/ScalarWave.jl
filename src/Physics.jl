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
        else # FIXME: Regression tests failing
	        operator[index] = (2.0/M)*chebw(i,N)*chebw(k,N)*chebd(k,l,N)*chebd(i,j,N)	
        end
	end
	return operator
end

function operatorNBC{T<:Int}(N::T, M::T)::Array{Float64, 4}
	operator = zeros(N+1, N+1, N+1, N+1)
    for index in CartesianRange(size(operator))    
        k = index.I[1]
        i = index.I[2]
        l = index.I[3]
        j = index.I[4]   
	    operator[index] = (2.0/M)*chebw(i,N)*chebw(k,N)*chebd(k,l,N)*chebd(i,j,N)	
	end
	return operator
end

function operatorNBCW{T<:Int}(N::T, M::T)::Array{Float64, 4}
	operator = zeros(N+1, N+1, N+1, N+1)
    for index in CartesianRange(size(operator))    
        k = index.I[1]
        i = index.I[2]
        l = index.I[3]
        j = index.I[4]   
	    operator[index] = (2.0/M)*chebd(k,l,N)*chebd(i,j,N)	
	end
	return operator
end

function initializeRHS{T<:Integer}(N::T, M::T, loc::Array{T,1}, fnbrow::Function, fnbcol::Function)::Array{Float64,2}	
	B = zeros(N+1, N+1)
    for index in CartesianRange(size(B))
        i = index.I[1]
        j = index.I[2]
		xp, yp = coordtrans(M, [chebx(i,N),chebx(j,N)], loc)
		if i == 1
			B[index] = fnbrow(xp,yp) 
        elseif j==1
            B[index] = fnbcol(xp,yp)
        else    
            continue
		end
	end
	return B
end
