#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2019
# Non-Linear Solver Routines
#--------------------------------------------------------------------

function linearOP(f, r, ϕ)::Array{Float64,2}
    return [vec(linearH(f, r, ϕ, :H1, :Δf)) vec(linearH(f, r, ϕ, :H1, :Δr)) vec(linearH(f, r, ϕ, :H1, :Δϕ));
            vec(linearH(f, r, ϕ, :H2, :Δf)) vec(linearH(f, r, ϕ, :H2, :Δr)) vec(linearH(f, r, ϕ, :H2, :Δϕ));
            vec(linearH(f, r, ϕ, :H3, :Δf)) vec(linearH(f, r, ϕ, :H3, :Δr)) vec(linearH(f, r, ϕ, :H3, :Δϕ))]
end

function boundaryOP(f, r, ϕ)::Array{Float64,2}
    return [vec(boundaryOP(f, r, ϕ, :Δf)) vec(boundaryOP(f, r, ϕ,  :0)) vec(boundaryOP(f, r, ϕ,  :0));
            vec(boundaryOP(f, r, ϕ,  :0)) vec(boundaryOP(f, r, ϕ, :Δr)) vec(boundaryOP(f, r, ϕ,  :0));
            vec(boundaryOP(f, r, ϕ,  :0)) vec(boundaryOP(f, r, ϕ,  :0)) vec(boundaryOP(f, r, ϕ, :Δϕ))]
end

function linearRHS(f, r, ϕ)::Array{Float64,1}
    return [vec(rhsH(:H1));
            vec(rhsH(:H2));
            vec(rhsH(:H3))]
end

function norm(f, r, ϕ)::Float64 
    return  norm([norm(linearH(f, r, ϕ, :H1, :Δf));
                  norm(linearH(f, r, ϕ, :H2, :Δf));
                  norm(linearH(f, r, ϕ, :H3, :Δf))])
end

function solve(S::T, A::Array{Float64,2}, b::Array{Float64,1}) where {T<:Cardinal{Tag}} where {Tag}
    x  = reshape(A \ b, (3, prod(order(S))))
    Δf = Field(S, x[1,:])
    Δr = Field(S, x[2,:])
    Δϕ = Field(S, x[3,:])
    return (Δf, Δr, Δϕ) 
end

function Newton(f0,  r0,  ϕ0, 
                fBC, rBC, ϕBC, maxiterations::Int, abstol::Float64)
    (f, r, ϕ) = (f0, r0, ϕ0)
    errornorm = norm(f, r, ϕ)
    for iteration in 1:maxiterations
        if errornorm < abstol
            return (f, r, phi)
        end
        # XXX: Check with Nathan on how to set BCs consistently
        (Δf, Δr, Δϕ) = solve(linearOP(f, r, ϕ) + boundaryOP(f, r, ϕ), linearRHS(f, r, ϕ))
        (f, r, ϕ) = (f + Δf, r + Δr, ϕ + Δϕ)
        setBCs(f, r, ϕ, fBC, rBC, ϕBC)
        errornorm = norm(f, r, ϕ)
   end
   @warn "Iterations did not converge. \n 
          Try with a different tolerance, or use a different initial guess"
   return 
end

