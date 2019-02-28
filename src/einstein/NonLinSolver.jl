#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2019
# Non-Linear Solver Routine
# #--------------------------------------------------------------------

function Newton(Svec, Fvec, Jvec, maxiterations::Int, abstol::Float64)
    for iteration in 1:maxiterations
        error = norm(Svec)
        if error < abstol
            println("The solver converged at iteration %i with residual %f", iteration, error)
            return Svec0 
        end
        ΔSvec = (Jvec + Bvec) \ Fvec
        Svec  = Svec + ΔSvec 
   end
   @warn "The solver failed to converge in $maxiterations iterations. Currest residual $error" 
end

