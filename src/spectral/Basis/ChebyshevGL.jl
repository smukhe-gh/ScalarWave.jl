#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 05-2019
# Chebyshev Polynomials: Gauss Lobatto grid 
# See Boyd F.8
#--------------------------------------------------------------------

function collocation(S::ChebyshevGaussLobatto{Tag, N, min, max}, i::Int) where {Tag, N, min, max}
    @assert max > min
    @assert i <= N+1
    return cospi(abs((i-1))/N)*(max - min)/2 + (max + min)/2
end

function derivative(S::ChebyshevGaussLobatto{Tag, N, min, max}, i::Int, j::Int) where {Tag, N, min, max}
    @assert i <= N+1
    @assert j <= N+1
	if i==j==1
		return (2N^2 + 1)/6
	elseif i==j==N+1
		return -(2N^2 + 1)/6
	elseif i==j
		return -chebx(j, N)/(2(1-chebx(j, N)^2))
	else
		ci = (i == 1 || i == N+1) ? 2 : 1
		cj = (j == 1 || j == N+1) ? 2 : 1
		s  = (i + j) % 2 != 0 ? -1 : 1
		return (ci/cj)*(s/(chebx(i,N)-chebx(j,N)))*(2/(max - min))
	end
end

function quadrature(S::ChebyshevGaussLobatto{Tag, N, min, max}, i::Int) where {Tag, N, min, max}
    # FIXME: Check order required for exact integration
	W = 0.0
	for j in 1:N+1
		w = (j == 1 ? 1 : (j-1)%2 == 0 ? 2/(1-(j-1)^2) : 0)
		l = (i == 1 || i == N+1 ? (1/N)*cospi((i-1)*(j-1)/N) : (2/N)*cospi((i-1)*(j-1)/N))
		W +=  w*l
	end
	return W*(max - min)/2

end
