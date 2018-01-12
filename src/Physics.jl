#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

function operator{T<:Integer}(N::T)	
	operator = eye((N+1)^2, (N+1)^2)
	for i in 1:N+1, j in 1:N+1, k in 1:N+1, l in 1:N+1
		if 	!(i==1 || k==1)
			I = (N+1)*(i-1) + k
			J = (N+1)*(j-1) + l
			for m in 1:N+1, n in 1:N+1
				operator[I, J] += w[i]*w[k]*(delta[i,m]*D[k,n]*delta[m,j]*D[n,l]) 
			end
		end
	end
	return operator
end

function initializeB{T<:Integer}(N::T, loc::Array{T,1}, M::T)	
	b = zeros((N+1)^2)
	bl(x,y) = sin(pi*x)
	br(x,y) = sin(pi*y)
	ve(x,y) = 0
	for i in 1:(N+1), j in 1:(N+1)
		I = (N+1)*(i-1) + j
		xp, yp = coortrans([x,y], loc, M)
		if i == 1
			b[I] = bl(xp, yp)
		elseif j == 1
			b[I] = br(xp, yp)
		else
			b[I] = w[i]*w[j]*ve(xp, yp)
		end
	end
	return b
end



