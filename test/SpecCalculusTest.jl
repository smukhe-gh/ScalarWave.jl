#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

@test chebx(3, 10) ≈ 0.8090169943749475
@test chebd(3, 1, 5) ≈ -0.723606797749979
@test chebw(9, 12) ≈ 0.2258075258075258

@test chebgrid(23)[7] == chebx(7,23)
@test chebweights(12)[8] == chebw(8,12)
@test vandermonde(10, chebgrid(14))[3,7] == cheb(6, chebx(3,14)) 

function consistency_chebx(N::Int)
	x  = chebx(5, N)
	Tx = cos(N*acos(x))
	return abs(Tx)
end
@test consistency_chebx(rand(10:20)) ≈ 1.0

function consistency_chebd(N::Int)
	DR = Float64[chebd(fld(N, 3), j, N) for j in 1:N+1]
	return sum(DR)
end
@test isapprox(consistency_chebd(rand(10:20)), 0.0; atol = 15)

function consistency_chebw(N::Int)
	w = Float64[chebw(i, N) for i in 1:N+1]
	return sum(w)
end
@test consistency_chebw(rand(2:20)) ≈ 2.0
