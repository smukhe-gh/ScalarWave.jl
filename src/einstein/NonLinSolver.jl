#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2019
# Non-Linear Solver Routine
# #--------------------------------------------------------------------

function Newton(space::Type{S}, Xvec::Array{Float64,1},
                Fvec::Function, Jvec::Function, Bvec::Function, 
                maxiterations::Int, abstol::Float64)::Array{Float64,1} where {S<:Space{Tag}} where {Tag}
    for iteration in 1:maxiterations
        error = norm(Fvec(space, Xvec))
        @show iteration, error
        if error < abstol
            println("The solver converged at iteration $iteration with residual $error")
            return Xvec 
        end
        ΔXvec = (Jvec(space, Xvec) + Bvec(space, Xvec)) \ -Fvec(space, Xvec)
        Xvec  = Xvec + ΔXvec 
   end
   @warn "The solver failed to converge in $maxiterations iterations." 
   return Xvec
end

