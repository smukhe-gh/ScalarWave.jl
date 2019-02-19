#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2019
# Non-Linear Solver Routines
#--------------------------------------------------------------------

function linearOP(f::Field{S}, r::Field{S}, ϕ::Field{S})::Array{Float64,2} where {S}
    return [vec(linearH(f, r, ϕ, :H1, :Δf)) vec(linearH(f, r, ϕ, :H1, :Δr)) vec(linearH(f, r, ϕ, :H1, :Δϕ));
            vec(linearH(f, r, ϕ, :H2, :Δf)) vec(linearH(f, r, ϕ, :H2, :Δr)) vec(linearH(f, r, ϕ, :H2, :Δϕ));
            vec(linearH(f, r, ϕ, :H3, :Δf)) vec(linearH(f, r, ϕ, :H3, :Δr)) vec(linearH(f, r, ϕ, :H3, :Δϕ))]
end

function boundaryOP(f::Field{S}, r::Field{S}, ϕ::Field{S})::Array{Float64,2} where {S}
    return [vec(boundary(f, r, ϕ, :Δf)) vec(boundary(f, r, ϕ,  :z)) vec(boundary(f, r, ϕ,  :z));
            vec(boundary(f, r, ϕ,  :z)) vec(boundary(f, r, ϕ, :Δr)) vec(boundary(f, r, ϕ,  :z));
            vec(boundary(f, r, ϕ,  :z)) vec(boundary(f, r, ϕ,  :z)) vec(boundary(f, r, ϕ, :Δϕ))]
end

function linearRHS(f::Field{S}, r::Field{S}, ϕ::Field{S})::Array{Float64,1} where {S}
    return [vec(rhsH(f, r, ϕ, :H1));
            vec(rhsH(f, r, ϕ, :H2));
            vec(rhsH(f, r, ϕ, :H3))]
end

function LinearAlgebra. norm(f::Field{S}, r::Field{S}, ϕ::Field{S})::Float64 where {S}
    return  norm([norm(H(f, r, ϕ, :H1));
                  norm(H(f, r, ϕ, :H2));
                  norm(H(f, r, ϕ, :H3))])
end

function LinearAlgebra. norm(f::Field{S})::Float64 where {S}
    return norm(vec(f));
end

function delta(f::Field{S}, r::Field{S}, ϕ::Field{S}) where {S}
    A  = linearOP(f, r, ϕ) .+ boundaryOP(f, r, ϕ)
    b  = linearRHS(f, r, ϕ)
    @show cond(A)
    x  = reshape(A \ b, (3, Int(length(b)/3)))
    Δf = Field(S, shape(S, x[1,:]))
    Δr = Field(S, shape(S, x[2,:]))
    Δϕ = Field(S, shape(S, x[3,:]))
    return (Δf, Δr, Δϕ) 
end

function Newton(f0,  r0,  ϕ0, 
                fBC, rBC, ϕBC, maxiterations::Int, abstol::Float64)
    (f, r, ϕ) = (f0, r0, ϕ0)
    errornorm = norm(f, r, ϕ)
    for iteration in 1:maxiterations
        @show iteration, errornorm
        if errornorm < abstol
            return (f, r, phi)
        end
        (Δf, Δr, Δϕ) = delta(f, r, ϕ)
        (f, r, ϕ) = (f + Δf, r + Δr, ϕ + Δϕ)
        (f, r, ϕ) = setBCs(f, r, ϕ, fBC, rBC, ϕBC)  # is this neccessary?
        errornorm = norm(f, r, ϕ)
   end
   @warn "Iterations did not converge. \n Try with a different tolerance, or use a different initial guess"
   return (f, r, ϕ)
end

