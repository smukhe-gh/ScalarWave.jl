#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

function operator{T<:Int}(N::T, M::T)::Array{Float64, 4}
	operator = zeros(N+1, N+1, N+1, N+1)
	#for k in 1:N+1, i in 1:N+1, l in 1:N+1, j in 1:N+1
    for index in CartesianRange(size(operator))    
        k = index.I[1]
        i = index.I[2]
        l = index.I[3]
        j = index.I[4]   
		
        if 	i==1 || k==1
			#operator[k, i, l, j] = delta(i,j)*delta(k,l)
		    operator[index] = delta(i,j)*delta(k,l)
        else
			#operator[k, i, l, j] = (2.0/M)*chebw(i,N)*chebw(k,N)*chebd(k,l,N)*chebd(i,j,N)
	        operator[index] = (2.0/M)*chebw(i,N)*chebw(k,N)*chebd(k,l,N)*chebd(i,j,N)	
        end
	end
	return operator
end

function setB{T<:Integer}(N::T, M::T, loc::Array{T,1})::Array{Float64,2}	
	B = zeros(N+1, N+1)

    bl(x,y) = 1.0
	br(x,y) = 1.0
	v(x,y)  = 0.0

	for i in 1:(N+1), j in 1:(N+1)
		xp, yp = coortrans([chebx(i),chebx(j)], loc, M)
		if i == 1
			B[j,i] = bl(xp, yp)
		elseif j == 1
			B[j,i] = br(xp, yp)
		else
			B[j,i] = chebx(i, N)*chebx(j,N)*v(xp, yp)
		end
	end
	return B
end



