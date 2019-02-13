#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2019
# Bad choice for a file name. This holds the expressions for
# linear and non-operators
# There are 4 equations: two hyperbolic in nature (H1, H2)
# and two elliptic in nature (E3, E4). 
# See Gundlach and Pullin 1997
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Compute residuals for non-linear equations
#--------------------------------------------------------------------

function H1(f::Field{S}, r::Field{S}, ϕ::Field{S})::Field{S} where {S<:Cardinal{Tag}} where {Tag}
    DV, DU = derivative(S)
    return DV*DU*r + ((DU*r)*(DV*r))/r - f/r
end

function H2(f::Field{S}, r::Field{S}, ϕ::Field{S})::Field{S} where {S<:Cardinal{Tag}} where {Tag}
    DV, DU = derivative(S)
    return DV*DU*log(f) + (2/r)*DV*DU*r + 2*(DU*ϕ)*(DV*ϕ)
end

function E1(f::Field{S}, r::Field{S}, ϕ::Field{S})::Field{S} where {S<:Cardinal{Tag}} where {Tag}
    DV, DU = derivative(S)
    return DU*DU*r - ((DU*f)/f)*DU*r + (DU*ϕ)^2
end

function E2(f::Field{S}, r::Field{S}, ϕ::Field{S})::Field{S} where {S<:Cardinal{Tag}} where {Tag}
    DV, DU = derivative(S)
    return DV*DV*r - ((DV*f)/f)*DV*r + (DV*ϕ)^2
end

#--------------------------------------------------------------------
# Compute linear operators 
#--------------------------------------------------------------------

function linearH1(f::Field{S}, r::Field{S}, 
                  var::Symbol)::ProductSpaceOperator{S} where {S<:Cardinal{Tag}} where {Tag}
    DV, DU = derivative(S)
    I      = identity(S)
    if var == :Δf
        return -1/r   
    else
        return DV*DU + ((DV*r)/r)*DU + ((DU*r)/r)*DV + (f/r^2)*I - ((DU*r)*(DV*r)/r^2)*I
    end
end

function linearH2(f::Field{S}, r::Field{S})::ProductSpaceOperator{S} where {S<:Cardinal{Tag}} where {Tag}
    DV, DU = derivative(S)
    I      = identity(S)
    if var == :Δf
        return (DU*DV/f) - ((DV*f)/f^2)*DU + ((DU*f)/f^2)*DV + 2*((DU*f)*(DV*f)/f^3)*I - 2*((DU*DV*r)/r^2)*I

    else
        return (2/r)*DV*DU + (2/r^2)*(DU*DV*r)
    end
end

function linearE1(f::Field{S}, r::Field{S})::Field{S} where {S<:Cardinal{Tag}} where {Tag}
    if var == :Δf
        return ((DU*r)/f)*DU + ((DU*f)*(DU*r)/f^2)*I
    else
        return DU*DU - ((DU*f)/f)*DU
    end
end

function linearE2(f::Field{S}, r::Field{S})::Field{S} where {S<:Cardinal{Tag}} where {Tag}
    if var == :Δf
        return ((DV*r)/f)*DV + ((DV*f)*(DV*r)/f^2)*I
    else
        return DV*DV - ((DV*f)/f)*DV
    end
end

function linearH2(f::Field{S}, r::Field{S})::Field{S} where {S<:Cardinal{Tag}} where {Tag}
    return DU*DV + ((DU*r)/r) + (DV*r)/r
end

#--------------------------------------------------------------------
# Now construct the whole operator and the Newton Iterator 
#--------------------------------------------------------------------

function lsolve(f::Field{S}, r::Field{S})
    # This operator will have a kernel
    # Where to impose these boundary conditions?
    linearoperator = [[linearH1(f, r, :Δf), linearH1(f, r, :Δr)],
                      [linearH2(f, r, :Δf), linearH2(f, r, :Δr)]]
    return solve(linearopeator, [Δf, Δr])
end

function initialdata(S::ProductSpace{GaussLobatto{Tag1, N1, max1, min1},
                                     GaussLobatto{Tag2, N2, max2, min2}}, mapf::Function, 
                                                                          mapr::Function)
    return (Field(S, mapf), Field(S, mapf))
end

# N.B. Define fmap, rmap, and ϕBC
f0, r0 = initialconditions(fmap, rmap)
ϕ0     = solve(linearH2(f0, r0), ϕBC)

function Newton(f0::Field{S}, r0::Field{S}, ϕ0::Field{S}, abstol::Float64, maxiter::Int)
    resE1 = E1(f0, r0, ϕ0)
    resE2 = E2(f0, r0, ϕ0)
    iter  = 0
    while norm([resE1, resE2]) > abstol $$ iter < maxiter
        [Δf, Δr] = lsolve(f, r)
        [f,   r] + [f, r] .+ [Δf, Δr]
        ϕ        = solve(linearH2(f, r), ϕBC)
        (resE1, resE2) = (E1(f, r, ϕ), E2(f, r, ϕ))
        iter     = iter + 1
    end
end
