#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2019
# Get the fucking coefficents! [Done!]
# Now check if they work in arbitrary spaces
#--------------------------------------------------------------------

function prefactor(i, N)
    (i == 1 || i == N+1) ? (return 1/2) : (return 1)
end

# construct the function
N = 3
c = rand(4) 
f = [sum(prefactor(m,N)*cheb(m-1, chebx(i, 3))*c[m] for m in 1:4) for i in 1:4]

# now get the coefficents back
r = (2/N).*[sum(prefactor(i,N)*cheb(m-1, chebx(i, 3))*f[i] for i in 1:4) for m in 1:4]
@test c ≈ r

# Now, try the same with spaces

# construct the function
function transform(α::Field{Chebyshev{Tag, N, max, min}})::Field{GaussLobatto{Tag, N, max, min}} where {Tag, N, max, min}
    u = Field(GaussLobatto{Tag, N, max, min})
    for gridindex in 1:N+1
        u.value[gridindex] = sum(prefactor(order+1, N)*cheb(order, chebx(gridindex, N))*α.value[order+1] for order in 0:N)
    end
    return u
end

# now get the coefficents back
function transform(u::Field{GaussLobatto{Tag, N, max, min}})::Field{Chebyshev{Tag, N, max, min}} where {Tag, N, max, min}
    α = Field(Chebyshev{Tag, N, max, min})
    for order in 0:N
        α.value[order+1] = (2/N)*sum(prefactor(gridindex, N)*cheb(order, chebx(gridindex, N))*u.value[gridindex] for gridindex in 1:N+1)
    end
    return α
end

S  = Chebyshev(U, 3)
cc = Field(S, c)
ff = transform(cc)
rr = transform(ff)

@test rr ≈ cc
