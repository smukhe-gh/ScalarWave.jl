#--------------------------------------------------------------------
# Solve a 1D non-linear equation.  
# u_xx = Exp[u]; u(+/-1) =  0
#--------------------------------------------------------------------

S = GaussLobatto(U, 36)
x = Field(S, x->0)
abstol  = 1e-10
maxiter = 40

#--------------------------------------------------------------------
# The non-linear and linearized equations
#--------------------------------------------------------------------

function F(u::Field{S})::Field{S} where {S}
    D = derivative(S)
    return D*D*u - exp(u)
end

function J(u::Field{S}, Δu::Union{Symbol, Field{S}})::Operator{S} where {S}
    D = derivative(S)
    return D*D*Δu - exp(u)*Δu
end


#--------------------------------------------------------------------
# Interface with the Newton solver 
#--------------------------------------------------------------------
#
function Fvec(space::S, uvec::Array{Float64,1})::Array{Float64,1} where {S}
    return vec(F(shape(space, uvec))) 
end

function Jvec(space::S, uvec::Array{Float64,1})::Array{Float64,2} where {S}
    return vec(J(shape(space, uvec), :u))
end

function Bvec(space::S, uvec::Array{Float64,1})::Array{Float64,2} where {S}
    return vec(boundary(space))
end

xsolved = Newton(S, vec(x), 
                 Fvec, Jvec, Bvec, 
                 maxiter, abstol) 
