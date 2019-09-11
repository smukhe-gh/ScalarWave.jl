#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 10-2019
# Test the performance of the residual evaluator
#--------------------------------------------------------------------
# [1] Compute operators correctly.
# [2] Compute product of operators and fields correctly.
# [3] Compute product of fields correctly. 
# [4] Can you simplify [1] and [2] into a single step? 


function F(a::Field{S}, r::Field{S}, ϕ::Field{S}, DU::Operator{S}, DV::Operator{S})::Field{S} where {S}
    F1 = r*(DU*(DV*ϕ)) + (DU*r)*(DV*ϕ) + (DV*r)*(DU*ϕ)
    return F1
end

function F12(a::Field{S}, r::Field{S}, ϕ::Field{S}, DU::Operator{S}, DV::Operator{S})::NTuple{3, Field{S}} where {S}
    F1 = r*(DU*(DV*ϕ)) + (DU*r)*(DV*ϕ) + (DV*r)*(DU*ϕ)
    F2 = r*(DU*(DV*r)) + (DU*r)*(DV*r) + (1/4)*(a^2)
    F3 = (1/a)*(DU*(DV*a)) - (1/a^2)*(DU*a)*(DV*a) + (1/r)*(DU*(DV*r)) + 4pi*(DU*ϕ)*(DV*ϕ)
    return (F1, F2, F3)
end

function unloop0(DU, I2, D1)
    A = zeros(size(DU.value)) 
    for index in CartesianIndices(A)
        i, j, k, l = index.I
        A[index]   = D1.value[i,k] * I2.value[j,l]
    end
    return A
end

function unloop1(DU, I2, D1, u)
    w = zeros(size(u.value))
    for index in CartesianIndices(w)
        i, j = index.I
        w[index] = sum(D1.value[i,k] * I2.value[j,l] * u.value[k,l] for k in 1:N, l in 1:N+1)
    end
    return w
end

function unloop2(DU, I2, D1, u, v)
    w = zeros(size(u.value))
    for index in CartesianIndices(w)
        i, j = index.I
        w[index] = sum(v.value[i,j] * D1.value[i,k] * I2.value[j,l] * u.value[k,l] for k in 1:N, l in 1:N+1)
    end
    return w
end

function unloop3(DU, DV, I1, I2, D1, D2)
    A = zeros(size(DU.value)) 
    for index in CartesianIndices(A)
        i, j, k, l = index.I
        A[index]   = sum(D1.value[i, m] * I2.value[j, n] * I1.value[m, k] * D2.value[n, l] for m in 1:N, n in 1:N+1)
    end
    return A
end

function unloop4(DU, DV, I1, I2, D1, D2, u)
    w = zeros(size(u.value)) 
    for index in CartesianIndices(w)
        i, j = index.I
        w[index]   = sum(u.value[k,l] * D1.value[i, m] * I2.value[j, n] * I1.value[m, k] * D2.value[n, l] for m in 1:N, n in 1:N+1, k in 1:N, l in 1:N+1)
    end
    return w
end

function unloop5(DU, DV, I1, I2, D1, D2, u, v)
    w = zeros(size(u.value)) 
    for index in CartesianIndices(w)
        i, j = index.I
        w[index]   = sum(v.value[i,j] * u.value[k,l] * D1.value[i, m] * I2.value[j, n] * I1.value[m, k] * D2.value[n, l] for m in 1:N, n in 1:N+1, k in 1:N, l in 1:N+1)
    end
    return w
end

function unloop6(DU, DV, I1, I2, D1, D2, u, v)  # <-- optimized
    w = zeros(size(u.value)) 
    for index in CartesianIndices(w)
        i, j = index.I
        w[index]   = sum(v.value[i,j] * u.value[k,l] * D1.value[i, k] * D2.value[j, l] for k in 1:N, l in 1:N+1)
    end
    return w
end

function unloop7(DU, DV, I1, I2, D1, D2, u, v)
    w = zeros(size(u.value)) 
    for index in CartesianIndices(w)
        i, j = index.I
        w[index] = sum(D1.value[i,k] * I2.value[j,l] * u.value[k,l] * I1.value[i,m] * D2.value[j,n] * v.value[m,n] for k in 1:N, l in 1:N+1, m in 1:N, n in 1:N+1)
    end
    return w
end

function unloop8(DU, DV, I1, I2, D1, D2, u, v) # <-- optimized
    w = zeros(size(u.value)) 
    for index in CartesianIndices(w)
        i, j = index.I
        w[index] = sum(D1.value[i,k] * u.value[k,j] * D2.value[j,n] * v.value[i,n] for k in 1:N, n in 1:N+1)
    end
    return w
end

function unloop9(a::Field{S}, r::Field{S}, ϕ::Field{S}, DU::Operator{S1}, DV::Operator{S2})::Field{S} where {S, S1, S2}
    F1 = zeros(size(a.value))
    for index in CartesianIndices(F1)
        i, j = index.I
        F1[index] = (sum(r.value[i,j] * ϕ.value[k,l] * DU.value[i, k] * DV.value[j, l] for k in 1:N, l in 1:N+1)
                    + sum(DU.value[i,k] * r.value[k,j] * DV.value[j,n] * ϕ.value[i,n] + DU.value[i,k] * ϕ.value[k,j] * DV.value[j,n] * r.value[i,n] for k in 1:N, n in 1:N+1))
    end
    return Field(a.space, F1)
end

function unloop10(a::Field{S}, r::Field{S}, ϕ::Field{S}, DU::Operator{S1}, DV::Operator{S2})::Field{S} where {S, S1, S2} # <--optimized
    F1 = zeros(size(a.value))
    @inbounds for index in CartesianIndices(F1)
        i, j = index.I
        F1[index] = (sum(r.value[i,j] * ϕ.value[k,n] * DU.value[i, k] * DV.value[j,n] 
                         + DU.value[i,k] * r.value[k,j] * DV.value[j,n] * ϕ.value[i,n] 
                         + DU.value[i,k] * ϕ.value[k,j] * DV.value[j,n] * r.value[i,n] for k in 1:N, n in 1:N+1))
    end
    return Field(a.space, F1)
end

function unloop11(a::Field{S}, r::Field{S}, ϕ::Field{S}, DU::Operator{S1}, DV::Operator{S2})::Field{S} where {S, S1, S2} # <--optimized with partial sums
    F1 = r.value.*(DU.value*(ϕ.value*transpose(DV.value))) + (DU.value*r.value).*transpose((DV.value*transpose(ϕ.value))) + (DU.value*ϕ.value).*transpose((DV.value*transpose(r.value)))
    return Field(a.space, F1)
end

function unloop12(a::Field{S}, r::Field{S}, ϕ::Field{S}, DU::Operator{S1}, DV::Operator{S2})::NTuple{3, Field{S}} where {S, S1, S2} # <--optimized with partial sums
    F1 = r.value.*(DU.value*(ϕ.value*transpose(DV.value))) + (DU.value*r.value).*transpose((DV.value*transpose(ϕ.value))) + (DU.value*ϕ.value).*transpose((DV.value*transpose(r.value)))
    F2 = r.value.*(DU.value*(r.value*transpose(DV.value))) + (DU.value*r.value).*transpose((DV.value*transpose(r.value))) + (1/4)*(a.value.*a.value)
    F3 = (1/a).value.*(DU.value*(a.value*transpose(DV.value))) + (1/a^2).value.*(DU.value*a.value).*transpose((DV.value*transpose(a.value))) + 4pi*(DU.value*ϕ.value).*transpose((DV.value*transpose(ϕ.value)))
    return (Field(a.space, F1), Field(a.space, F2), Field(a.space, F3))
end
#----------------------------------------------------------
# Test loop expansions 
#----------------------------------------------------------

using Profile

struct U end
struct V end

N  = 10
PS = ProductSpace(ChebyshevGL{U, N, Float64}(-8, -6), ChebyshevGL{V, N+1, Float64}( 3,  5))
DU,DV = derivative(PS)
u  = Field(PS, (u,v)->u)
v  = Field(PS, (u,v)->v)
x  = Field(PS, (u,v)->u*v)
y  = Field(PS, (u,v)->u+v)

D1 = derivative(PS.S1)
D2 = derivative(PS.S2) 
I1 = identity(PS.S1)
I2 = identity(PS.S2) 

@test unloop0(DU, I2, D1) == DU.value
@test unloop1(DU, I2, D1, x) ≈ v.value
@test unloop2(DU, I2, D1, x, y) ≈ (y*v).value
@test unloop3(DU, DV, I1, I2, D1, D2) == (DU*DV).value
@test unloop4(DU, DV, I1, I2, D1, D2, x) ≈ (DU*(DV*x)).value
@test unloop5(DU, DV, I1, I2, D1, D2, x, v) ≈ v.value
@test unloop6(DU, DV, I1, I2, D1, D2, x, v) ≈ v.value
@test unloop7(DU, DV, I1, I2, D1, D2, x, x) ≈ (u*v).value
@test unloop8(DU, DV, I1, I2, D1, D2, x, x) ≈ (u*v).value

a = Field(PS, (u,v)->exp(u*v))
r = Field(PS, (u,v)->v-u)
ϕ = Field(PS, (u,v)->sin(u)*sin(v))

@test F(a, r, ϕ, DU, DV) ≈ unloop9(a, r, ϕ, D1, D2)
@test L2(F(a, r, ϕ, DU, DV) - unloop9(a, r, ϕ, D1, D2)) < 1e-12
@test L2(F(a, r, ϕ, DU, DV) - unloop10(a, r, ϕ, D1, D2)) < 1e-12
@test L2(F(a, r, ϕ, DU, DV) - unloop11(a, r, ϕ, D1, D2)) < 1e-12
# @test maximum.(F12(a, r, ϕ, DU, DV) .- unloop12(a, r, ϕ, D1, D2)) .<  1e-12
# exit()

#----------------------------------------------------------
# Now compare performance
#----------------------------------------------------------

N  = 100
PS = ProductSpace(ChebyshevGL{U, N, Float64}(-8, -6), ChebyshevGL{V, N+1, Float64}( 3,  5))
DU,DV = derivative(PS)
D1 = derivative(PS.S1)
D2 = derivative(PS.S2) 
I1 = identity(PS.S1)
I2 = identity(PS.S2) 
a  = Field(PS, (u,v)->exp(u*v))
r  = Field(PS, (u,v)->v-u)
ϕ  = Field(PS, (u,v)->sin(u)*sin(v))

@time F(a, r, ϕ, DU, DV) 
@time unloop11(a, r, ϕ, D1, D2)

println("\n==> F")
@timev F(a, r, ϕ, DU, DV) 

println("\n==> loopF")
@timev unloop11(a, r, ϕ, D1, D2)


