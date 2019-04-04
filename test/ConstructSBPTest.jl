#---------------------------------------------------------------
# Construct SBP operators for Gauss-Lobatto points 
# See <https://doi.org/10.1016/j.jcp.2014.01.038>
# Soham M 03/2019
#---------------------------------------------------------------
using LinearAlgebra

function constructH(w::Array{T,1})::Array{T,2} where {T}
    return diagm(0=>w)
end

function constructX(x)
    P = length(x) - 1
    X = zeros(P+1, P+1)
    for j in 0:P
        X[:,j+1] = x.^j
    end
    return X
end

function constructR(x,w)
    P = length(x) - 1
    R = zeros(P+1, P+1)
    H = constructH(w)
    Ẽ = constructẼ(x)
    for j in 0:P
        R[:,j+1] = j*H*(x.^(j-1)) - (1/2)*Ẽ*(x.^j)
    end
    return R
end

function constructΘ(x, w)
    Ẽ  = constructẼ(x)
    ΘS = Ẽ/2
    X  = constructX(x)
    R  = constructR(x,w)
    ΘA = R*inv(X)
    return (ΘS + ΘA)
end

function constructD(x, w)
    H = constructH(w)
    Θ = constructΘ(x,w)
    return inv(H)*Θ
end

function testSBP(x, w)
    H = constructH(w)
    Θ = constructΘ(x, w)  
    D = constructD(x, w)
    Ẽ = constructẼ(x)
    P = length(x) - 1
    @testset "SBP" begin
        for j in 1:P
            @test D*(x.^j) ≈ j*(x.^(j-1))
        end
        @test D == inv(H)*Θ
        @test Θ + transpose(Θ) ≈ Ẽ
    end
end

function constructẼ(x)
    P  = length(x) - 1
    ta = [basis(j, P, chebx(0,P)) for j in 0:P]
    tb = [basis(j, P, chebx(P,P)) for j in 0:P]
    e  = zeros(P)
    e0 = [1, e...]
    en = [e..., 1]
    T  = e0*ta' + en*tb'
    E  = zeros(P+1, P+1)
    E[1,1] = 1
    E[end, end] = -1
    Ẽ = T'*E*T
    return Ẽ
end

# Specify collocation points, quadrature weights and 
# associated basis functions to compute the SBP operator
# xx = [-1.0, -0.2898979485566356, 0.6898979485566357]
# w  = [0.2222222222222222, 1.0249716523768433, 0.7528061254009345]

# function basis(j,x)
    # S = 1.0
    # for m in 0:2
        # if m == j
            # S *= 1
        # else
            # S *= (x - xx[m+1])/(xx[j+1] - xx[m+1])
        # end
    # end
    # return S
# end

function chebx(i, N)
    return cospi(i/N)
end

function basis(j::Int, P::Int,  x::Float64)::Float64 
    S = 0
    for m in 0:P
        pm = (m == 0 || m == P) ? 2 : 1
        S += (1/pm)*cheb(m,x)*cheb(m, chebx(j, P))
    end
    pj = (j == 0 || j == P) ? 2 : 1
    return (2/(P*pj))*S
end

P  = 3
xx = [chebx(i, P) for i in 0:P]  
w  = [chebw(i, P) for i in 1:P+1]
# FIXME: Test failing for Chebyshev collocation points
testSBP(xx,w)
