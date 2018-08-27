#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
# Rationals Test
# Choose uniform collocation points in a cardinal basis (polynomial)
# and compute the derivative and the integral
#--------------------------------------------------------------------

function collocation{T<:Int}(i::T, N::T)::Rational
    return -(-1 + (2*(i-1)//N))
end

function derivative{T<:Int}(i::T, j::T, N::T)::Rational
    x = Rational[collocation(k, N) for k in 1:N+1]
    if i == j
        return sum((k==j ? 0 : (x[j] - x[k])^(-1)) for k in 1:N+1)
    else
        ai = prod((k==i ? 1 : (x[i] - x[k])) for k in 1:N+1)
        aj = prod((k==j ? 1 : (x[j] - x[k])) for k in 1:N+1)
        return ai/(aj*(x[i] - x[j]))
    end
end

N    = 5
grid = Rational[collocation(k, N) for k in 1:N+1]
D    = Rational[derivative(i, j, N) for i in 1:N+1, j in 1:N+1]
D2   = similar(D)

for i in 1:N+1, j in 1:N+1
    D2[i,j] = sum(D[i,k]*D[k,j] for k in 1:N+1)
end

#--------------------------------------------------------------------
# Test simple derivative
#--------------------------------------------------------------------
ϕ    = grid.^2
ψ    = similar(ϕ)

for i in 1:N+1
    ψ[i] = sum(D[i,k]*ϕ[k] for k in 1:N+1)
end

@test ψ == 2.*grid

#--------------------------------------------------------------------
# Solve a 1D elliptic PDE 
#--------------------------------------------------------------------
x = grid
B = rationalize.(zeros(N+1, N+1))
B[1,1] = B[end, end] = 1//1
b = rationalize.(zeros(N+1))
b[1]   =  1//20
b[end] = -1//20
u = (D2 + B) \ (x.^3 + b)
@test u == x.^5/20

#--------------------------------------------------------------------
# Solve a 2D elliptic PDE 
#--------------------------------------------------------------------

I   = rationalize.(eye(N+1))
Dx  = rationalize.(zeros((N+1), (N+1), (N+1), (N+1)))
Dy  = similar(Dx)
Φ   = Rational[xx^3 + yy^2 for xx in x, yy in x] 
dΦx = similar(Φ)
dΦy = similar(Φ)

for index in CartesianRange(size(Dx))
    i,j,k,l   = index.I
    Dx[index] = D[j,l]*I[i,k] 
    Dy[index] = I[j,l]*D[i,k] 
end

for index in CartesianRange(size(Φ))
    i, ii      = index.I
    dΦx[index] = sum(Dy[i, ii, j, jj]*Φ[j, jj] for j in 1:N+1, jj in 1:N+1)
    dΦy[index] = sum(Dx[i, ii, j, jj]*Φ[j, jj] for j in 1:N+1, jj in 1:N+1)
end

@test dΦy == Rational[2*yy for xx in x, yy in x]
@test dΦx == Rational[3*xx^2 for xx in x, yy in x]
