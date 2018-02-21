#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

function consistency_chebx(N::Int)
	x  = chebx(5, N)
	Tx = cos(N*acos(x))
	return abs(Tx)
end

function consistency_chebd(N::Int)
	DR = Float64[chebd(fld(N, 3), j, N) for j in 1:N+1]
	return sum(DR)
end

function consistency_chebw(N::Int)
	w = Float64[chebw(i, N) for i in 1:N+1]
	return sum(w)
end

function check_chebd(N::Int)::Float64
    x     = chebgrid(N)
    fx    = x.^5
    edfx  = 5.*x.^4
    chebD = Float64[chebd(i,j,N) for i in 1:N+1, j in 1:N+1]
    ndfx  = chebD*fx
    return maximum(abs.(ndfx - edfx))
end

@test chebx(3, 10) ≈ 0.8090169943749475
@test chebw(9, 12) ≈ 0.2258075258075258
@test chebgrid(23)[7] == chebx(7,23)
@test chebweights(12)[8] == chebw(8,12)
@test vandermonde(10, chebgrid(14))[3,7] == cheb(6, chebx(3,14)) 
@test consistency_chebx(rand(10:20)) ≈ 1.0
@test isapprox(consistency_chebd(rand(10:20)), 0.0; atol = 15)
@test consistency_chebw(rand(2:20)) ≈ 2.0
@test check_chebd(10) < 1e-13
@test_broken check_chebd(20) < 1e-13

