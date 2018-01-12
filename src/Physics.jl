#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

function operator{T<:Integer}(N::T)	
	# operator = eye((N+1)^2, (N+1)^2)
	operator = eye(N+1, N+1, N+1, N+1)
	# Take care of the order!
	for i in 1:N+1, j in 1:N+1, k in 1:N+1, l in 1:N+1
		if 	!(i==1 || k==1)
			# I = (N+1)*(i-1) + k
			# J = (N+1)*(j-1) + l
			# for m in 1:N+1, n in 1:N+1
			# 	# TODO: Replace with appropriate functions for D and w
			# 	# operator[I, J] += w[i]*w[k]*(delta[i,m]*D[k,n]*delta[m,j]*D[n,l] + 
			# 	# 							 D[i,m]*delta[k,n]*delta[m,j]*D[n,k]) 
			# 	operator[k, i, l, j] += w[i]*w[k]*(delta[i,m]*D[k,n]*D[m,j]*delta[n,l] + 
			# 								       D[i,m]*delta[k,n]*delta[m,j]*D[n,l]) 
			# end
			operator[k, i, l, j] += 2*w[i]*w[k]*D[k,l]*D[i,j] 
		end
	end
	operator = reshape(operator, ((N+1)^2, (N+1)^2))
	return operator
end

function initializeB{T<:Integer}(N::T, loc::Array{T,1}, M::T)	
	# b = zeros((N+1)^2)
	b = zeros(N+1, N+1)
	bl(x,y) = sin(pi*x)
	br(x,y) = sin(pi*y)
	ve(x,y) = 0
	for i in 1:(N+1), j in 1:(N+1)
		# I = (N+1)*(i-1) + j
		xp, yp = coortrans([x,y], loc, M)
		if i == 1
			b[j,i] = bl(xp, yp)
		elseif j == 1
			b[j,i] = br(xp, yp)
		else
			b[j,i] = w[i]*w[j]*ve(xp, yp)
		end
	end
	b = reshape(b, ((N+1)^2,))
	return b
end



