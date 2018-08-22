#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
# All operations are done with rationals to check 
# if the results are exact
#--------------------------------------------------------------------

function chebxRational{T<:Int}(i::T, N::T)::Rational
    return rationalize(BigInt, cospi((i-1)/N), 1e-15)
end

function chebdRational{T<:Int}(i::T, j::T, N::T)::Rational
	if i==j==1
		return (2N^2 + 1)//6
	elseif i==j==N+1
		return -(2N^2 + 1)//6
	elseif i==j
		return -chebxRational(j, N)//(2(1-chebxRational(j, N)^2))
	else
		ci = (i == 1 || i == N+1) ? 2 : 1
		cj = (j == 1 || j == N+1) ? 2 : 1
		s  = (i + j) % 2 != 0 ? -1 : 1
		return (ci//cj)*(s//(chebxRational(i,N)-chebxRational(j,N)))
	end
end

function chebwRational{T<:Int}(i::T, N::T)::Rational
	W = 0//1
	for j in 1:N+1
		w = (j == 1 ? 1 : (j-1)%2 == 0 ? 2//(1-(j-1)^2): 0)
        l = (i == 1 || i == N+1 ? (1//N)*rationalize(BigInt, cospi((i-1)*(j-1)/N), 
                 1e-15) : (2//N)*rationalize(BigInt, cospi((i-1)*(j-1)/N), 1e-15))
		W +=  w*l
	end
	return W
end

P = 10
x = Rational[chebxRational(i, P) for i in 1:P+1]
w = Rational[chebwRational(i, P) for i in 1:P+1]
D = Rational[chebdRational(i, j, P) for i in 1:P+1, j in 1:P+1]

@test_broken sum(w) == 2//1
@test Float64(sum(w)) == 2.0

L = Array{Rational,2}(P+1, P+1)
for i in 1:P+1, j in 1:P+1
    L[i,j] = sum(D[i,k]*D[k,j] for k in 1:P+1)
end

ϕ = Rational[3*x^2 for x in x] 
u = Rational[6//1 for x in x] 

@test typeof(L) == Array{Rational,2}
@test typeof(ϕ) == Array{Rational,1}

ψ = Array{Rational,1}(P+1)
for i in 1:P+1
    ψ[i] = sum(L[i,k]*ϕ[k] for k in 1:P+1)
end

@test typeof(ψ) == Array{Rational,1}
@test_broken ψ == u 
@test Float64.(ψ) ≈ Float64.(u) 
