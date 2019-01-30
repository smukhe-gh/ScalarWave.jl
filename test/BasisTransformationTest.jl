#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2019
# Basis transformation Test
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# test 1D spaces
#--------------------------------------------------------------------

f = Field(GaussLobatto{U, 10, 2, -3}, rand(11))
fbar = basistransform(basistransform(f))
@test f ≈ fbar

#--------------------------------------------------------------------
# Temporary tests
#--------------------------------------------------------------------
# Create a random coefficent array
M = 4
N = 9
c =  rand(M+1, N+1)

function transform(c)
    f = zeros(M+1, N+1)
    for i in 1:M+1, j in 1:N+1
        f[i,j] = sum(prefactor(m+1,n+1, M, N)*c[m+1,n+1]*cheb(m, chebx(i, M))*cheb(n, chebx(j,N)) for m in 0:M, n in 0:N) 
    end
    return f
end

function reversetransform(u)
    c = zeros(M+1, N+1)
    for m in 0:M, n in 0:N
        c[m+1,n+1] = (4/(M*N))*sum(prefactor(i,j, M, N)*u[i,j]*cheb(m, chebx(i, M))*cheb(n, chebx(j,N)) for i in 1:M+1, j in 1:N+1) 
    end
    return c
end

# compute f from c
f = transform(c)

# now get back the coefficent vector from f
cc = reversetransform(f)

# now test if they are equal
@test c ≈ cc

#--------------------------------------------------------------------
# test 2D spaces
#--------------------------------------------------------------------

S = ProductSpace{GaussLobatto{V, 9, 1, -1}, 
                 GaussLobatto{U, 4, 1, -1}}

uu = Field(S, (U,V)->U)
vv = Field(S, (U,V)->V)
ff = Field(S, (U,V)->U^2 + V^3)

cc  = basistransform(ff)
ffr = basistransform(cc) 
ccr = basistransform(ffr) 
@test ffr.value ≈ ff.value
@test ccr.value ≈ cc.value
