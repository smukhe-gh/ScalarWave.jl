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

# evaluate a field at a single point
function evaluate(u::Field{GaussLobatto{Tag, N, max, min}}, x::Number)::Number where {Tag, N, max, min}
    α = transform(u)
    xmap = (x - (max + min)/2)*(2/(min - max))
    # XXX: Why do we need a minus sign?   
    return - sum(prefactor(order+1, N)*cheb(order, xmap)*α.value[order+1] for order in 0:N)
end

# check for domain from -1 to 1
c  = rand(4)
S  = Chebyshev(U, 3)
cc = Field(S, c)
ff = transform(cc)
rr = transform(ff)
@test rr ≈ cc

# Now check what happens when the domain is no longer -1 to 1
S  = Chebyshev(U, 3, 5, -3)
cc = Field(S, c)
ff = transform(cc)
rr = transform(ff)
@test rr ≈ cc

# now test the function evaluation at a point
f = Field(GaussLobatto{U, 10, 1, -1}, x->x^3)
@show evaluate(f, -0.8)

f = Field(GaussLobatto{U, 10, 5, -3}, x->x^3)
@show evaluate(f, -0.8)

f = Field(GaussLobatto{U, 10, 100, -30}, x->x^3)
@show evaluate(f, -0.8)

@show f(-0.8)
