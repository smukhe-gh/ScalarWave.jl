#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

function testchebx(N::Int)
	x  = chebx(5, N)
	Tx = cos(N*acos(x))
	return abs(Tx)
end

function testchebd(N::Int)
	DR = Float64[chebd(fld(N, 3), j, N) for j in 1:N+1]
	return sum(DR)
end

function testchebdpoly(N::Int)::Float64
    x     = chebgrid(N)
    fx    = x.^5
    edfx  = 5.*x.^4
    chebD = Float64[chebd(i,j,N) for i in 1:N+1, j in 1:N+1]
    ndfx  = chebD*fx
    return maximum(abs.(ndfx - edfx))
end

function testchebw(N::Int)
	w = Float64[chebw(i, N) for i in 1:N+1]
	return sum(w)
end

@test chebx(3, 10) ≈ 0.8090169943749475
@test chebw(9, 12) ≈ 0.2258075258075258
@test vandermonde(10)[2,5] == cheb(4, chebx(2,10))
@test pseudovandermonde(10, chebgrid(14))[3,7] == cheb(6, chebx(3,14)) 
@test chebgrid(23)[7] == chebx(7,23)
@test chebweights(12)[8] == chebw(8,12)

@test testchebx(rand(10:20)) ≈ 1.0
@test testchebd(rand(10:20)) < 1e-14
@test testchebw(rand(2:20)) ≈ 2.0
@test chebgrid(4, 2, 1) == (chebgrid(4) + 1)/2
@test chebgrid(4, 2, 2) == (chebgrid(4) - 1)/2
@test chebgrid(2, 2, 1)[end] == chebgrid(2, 2, 2)[1]
@test testchebdpoly(10) < 1e-13
