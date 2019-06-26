#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 05-2019
# Chebyshev Polynomials: Gauss Lobatto grid 
# See Boyd F.8
# Check if all type conversions are consistent with working precision
#--------------------------------------------------------------------

export collocation, derivative, integral 

function naturalcollocation(S::ChebyshevGL{Tag, N, T}, i::Int)::T where {Tag, N, T}
    @assert i <= N
    return  cospi(T(i-1)/T(N-1))
end


function naturalderivative(S::ChebyshevGL{Tag, N, T}, i::Int, j::Int)::T where {Tag, N, T}
    @assert i <= N && j <= N
	if i==j==1
        return  (2T(N-1)^2 + 1)/T(6)
	elseif i==j==N
        return -(2T(N-1)^2 + 1)/T(6)
	elseif i==j
        return - naturalcollocation(S,j)/(2(1-naturalcollocation(S,j)^2))
	else
		ci = (i == 1 || i == N) ?  2 : 1
		cj = (j == 1 || j == N) ?  2 : 1
		s  = (i + j) % 2 != 0   ? -1 : 1
        return (T(ci)/T(cj))*(T(s)/(naturalcollocation(S, i) - naturalcollocation(S,j)))
	end
end

function naturalintegral(S::ChebyshevGL{Tag, N, T}, i::Int)::T where {Tag, N, T}
    @assert i <= N
    W = T(0)
	for j in 1:N
        w = (j == 1 ? 1 : (j-1)%2 == 0 ? T(2)/T((1-(j-1)^2)) : 0)
        l = (i == 1 || i == N ? (T(1)/T(N-1))*cospi(T(i-1)*T(j-1)/T(N-1)) : (T(2)/T(N-1))*cospi(T(i-1)*T(j-1)/T(N-1)))
        W = W + w*l
	end
	return W
end

function collocation(space::S, i::Int)::T where {S <: Cardinal{Tag, N, T}} where {Tag, N, T}
    return naturalcollocation(space,i)*T(space.max - space.min)/2 + T(space.max + space.min)/2 
end

function derivative(space::S, i::Int, j::Int)::T where {S <: Cardinal{Tag, N, T}} where {Tag, N, T}
    return naturalderivative(space,i,j)*(2/T(space.max - space.min))
end

function integral(space::S, i::Int)::T where {S <: Cardinal{Tag, N, T}} where {Tag, N, T}
    return naturalintegral(space,i)*T(space.max - space.min)/T(2)
end
