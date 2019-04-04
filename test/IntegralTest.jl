#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 04-2018
# Test integration on a patch
#--------------------------------------------------------------------

function Base. *(A::IntegrationOperator{S}, B::Operator{S})::Operator{S} where {S}
    AP = Operator(S, A.value)  
    return AP*B
end

function Base. transpose(A::Operator{S})::Operator{S} where {S}
    B = transpose(A.value) .+ 0
    return Operator(S, B)
end

using LinearAlgebra

function Base. show(A::Operator{S}) where {S}
    display(A.value)
    println()
end

function Base. show(A::IntegrationOperator{S}) where {S}
    display(A.value)
    println()
end

function Base. inv(A::IntegrationOperator{S})::IntegrationOperator{S} where {S}
    return IntegrationOperator(S, inv(A.value))
end

function Base. transpose(u::Field{S}) where {S}
    return Field(S, transpose(u.value) .+ 0)
end

function dot(u::Field{S}, v::Field{S})::Float64 where {S}
    return sum((u*v).value) 
end

for P in 1:4
    S = GaussLobatto(U, P)
    D = derivative(S)
    # F = Field(S, x->x^4)
    # W = integral(S)  
    # @show P, W*F - 2/5
    show(D)
end

exit()

P = 4
S = GaussLobatto(U, P)
W = integral(S)
D = derivative(S)

show(W)
exit()
# We test
# <u,dw> + <w, du> = uLwL - uRwR 
# u'*W*D*w + w'*W*D*u = u'*W*B*w
# u'*W*D*w + u'(W*D)'*w = u'*W*B*w
B = inv(W)*(W*D + transpose(W*D))

# Test if B actually does project out the boundary
u = Field(S, x->x^(P-2))
v = Field(S, x->x^(P-1))

# @show dot(u, W*D*v) + dot(v, W*D*u)
# @show dot(u, W*D*v) + dot(u, transpose(W*D)*v)
# @show dot(u, W*B*v)
# show(D)
# show(B)


function quadratureweights()
