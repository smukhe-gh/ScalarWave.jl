#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 07-2019
# Write a non-linear solver from scratch for Einstein's field
# equations
#--------------------------------------------------------------------

function Base. log(u::Field{S})::Field{S} where {S}
    return Field(u.space, log.(u.value))
end

function J(f::Field{S}, r::Field{S}, ϕ::Field{S})::Tuple{Operator{S}} where {S}
    Jff = (- (r^2/(2*f^2))*(DU*DV) + (r^2/(2*f^3))*(DV*f)*DU + (r^2/(2*f^3))*(DU*f)*DV 
           - ((3*r^2)/(2*f^4))*(DU*f)*(DV*f)*I +  (r^2/f^3)*(DU*DV*f)*I + (r/f^2)*(DU*DV*r)*I)
    Jfr = - (r/f)*(DU*DV) + (r/f^3)*(DU*f)*(DV*f)*I - (r/f^2)*(DU*DV*f)*I - (1/f)*(DU*DV*r)*I
    Jfϕ = 0*I

    Jrf = (2/f^2)*I
    Jrr = (1/r)*(DU*DV) + (2/r^2)*(DV*r)*DU + (2/r^2)*DV - (4/r^3)*f*I - (4/r^3)*(DU*r)*(DV*r)*I - (2/r^2)*(DU*DV*r)*I
    Jrϕ = 0*I 

    Jϕf = 0*I 
    Jϕr = (1/r)*(DV*ϕ)*DU + (1/r)*(DV*ϕ)*DV - (1/r^2)*(DV*ϕ)*I 
    Jϕϕ = DU*DV + (1/r)*(DU*r)*DV + (1/r)*(DV*r)*DV
end

function *(op::Operator{S, T}, f::Field{S, Dual{T}})::Field{S, Dual{T}} where {S, T}
    return Dual(op * regular(f), [op * dual(f)[i] for i in 1:length(dual(f))])
end 

function JF(f::Field{S, Dual{T}})::Field{S, Dual{T}} where {S, T}
    resf = DU*DV*log(f)
    return resf
end

function F(f::Field{S, T})::Field{S, T} where {S, T}
    resf = DU*DV*log(f)
    return resf
end

function F(f::Field{S}, r::Field{S}, ϕ::Field{S})::Tuple{Field{S}} where {S}
    resf = DU*DV*log(f) + (2/r)*(DU*DV*r) + 2*(DU*ϕ)*(DV*ϕ) 
    resr = DU*DV*r +  (1/r)*(DU*r)*(DV*r) - (f/r)
    resϕ = DU*DV*ϕ + (1/r)*(DU*r)*(DV*ϕ) + (1/r)*(DV*r)*(DU*ϕ)
    return (bndf ⊕ resf, bndr ⊕ resr, bndϕ ⊕ resϕ)
end

function NewtonStep(u::Field{S}, J::Operator{S}, F::Field{S})::Field{S} where {S}
    Δ = solve(J, -F)
    return u + Δ
end

function NewtonIterator(u::Field{S}, J::Function, F::Function; maxiter=1000, abstol=1e-14, reltol=1e-14) where {S}
    println("Initial Newton residual = ", L2(F(u)))
    for iteration in 1:maxiter
        u = NewtonStep(u, J(u), F(u))
        println(iteration, "\t", L2(F(u)))
        L2(F(u)) < abstol ? (return u) : continue
    end
    return u 
end 
    
struct M end
S = ChebyshevGL{M, 15, BigFloat}(-1, 1)
D = derivative(S)
I = identity(S)
B = incomingboundary(S) ⊕ outgoingboundary(S)

s = Field(S, x->x)
b = Field(S, x->0)

u0 = (I-B)*s + B*b 
uf = NewtonIterator(u0, J, F, abstol=1e-20)

# using PyPlot
# @show uf.value[8]
# plot(uf)
# show()
