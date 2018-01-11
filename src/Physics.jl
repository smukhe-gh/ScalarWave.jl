#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

# FIXME: Shift to 'using' keyword for DArray/SharedArray to work.
include("SpecCalculus.jl") 
include("Utilities.jl")

module Physics

	# TODO: use CartesianIndices for this
	gidx(N, i0, i1) = i0-1 + (N+1) * (i1-1) + 1



	# Construct the entire operator element-wise, possibly using DArray/SharedArray.
	function operator(N::Int64)
		operator = zeros((N+1)^2, (N+1)^2)
		# TODO: or op = zeros(N+1, N+1, N+1, N+1)?
		#          M = reshape(op, ((N+1)^2, (N+1)^2))
		#          B = reshape(b, ((N+1)^2,))
		#          U = M \ B
		#          u = reshape(U, (N+1, N+1))
		#          op[i0,i1,j0,j1]
		#          M[i,j]   i=gidx(i0,i1), j=gidx(j0,j1)
		for i0 in 1:N+1, i1 in 1:N+1, j0 in 1:N+1, j1 in 1:N+1
		# for i in 1:(N+1)^2, j in 1:(N+1)^2	# Use Cartesian Indices
			i = gidx(N, i0, i1)
			j = gidx(N, j0, j1)
			if i0==i1==j0==j1==1
				println(SpecCalculus.chebw(i0, N)*SpecCalculus.chebw(j0, N))
			end
			if i0==1 || i1==1	# sets Dirichlet BCs [row stacked]
				operator[i,j] = Utilities.imat(i0, j0) * Utilities.imat(i1, j1)
			else # implement W[(I⊗Du)(Dv⊗I) + (Dv⊗I)(I⊗Du)] [FIXME]
				for k0 in 1:N+1, k1 in 1:N+1
					k = gidx(N, k0, k1)
					# XXX: You can eliminate storage completely by having a definition of op[i,j]
					operator[i,j] += SpecCalculus.chebw(i0, N)*SpecCalculus.chebw(i1, N)*
							(Utilities.imat(i0, k0)*SpecCalculus.chebd(k1, j1, N) + 
							SpecCalculus.chebd(i0, k0, N)*Utilities.imat(k1, j1))
				end
			end
		end
		return operator
	end

	function initialize(b::Array{Float64,2}, patchindex::Array{Int64,1}, M::Int64)
		
		bfuncl(x,y) = sin(pi*x)
		bfuncr(x,y) = sin(pi*y)
		vfunc(x,y)  = 0

		for i in 1:(N+1), j in 1:(N+1)	# User Cartesian Indices
			if i == 1
				b[i,j] = bfuncl(affine(chebx(i), M, d), affine(chebx(j), M, d))
			elseif j == 1
				b[i,j] = bfuncr(affine(chebx(i), M, d), affine(chebx(j), M, d))
			else
				b[i,j] = SpecCalculus.chebw(i, N)*SpecCalculus.chebw(j, N)*
						 vfunc(affine(chebx(i), patchindex[1], M), affine(chebx(j), patchindex[2], M))
			end
		end
		return b
	end

end 	# end module

