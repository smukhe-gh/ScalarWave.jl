#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 11-2018
# Compute the Vaidya metric for various mass functions
#   -- Use the analytic solutions to test
#   -- Compute scattering on the semi-analytic solutions
#--------------------------------------------------------------------


# The Vaidya metric for exponential mass function [Waugh & Lake, 1986]
function r_of_m_exp(u, v, α, β)
    c    = 1
    m(v) = (1/β)*(α*exp((β/2)*c*v) + 1)
    P(u) = - c*u
    x(r, v) = sqrt(r^2 - (4r/β) + (4*m(v)/β))
    f(r) = β*x(r, v) + 2*log(r - (2/β) + x(r,v)) - (β*c*v)/2 - P(u)
    r    = find_zero(f, 2) 
end


function f_of_m_exp(u, v, α, β)
    c    = 1
    m(v) = (1/β)*(α*exp((β/2)*c*v) + 1)
    P(u) = - c*u
    r    = r_of_m_exp(u, v, α, β)
    x    = sqrt(r^2 - (4r/β) + (4*m(v)/β)) 
    return (c^2*x)/(β*r)
end

# Mass function for a radial collapse of a null fluid [Saa and Girrotto, 2004]
# Use BigFloats for accuracy 
# TODO: Experiment with the ODE solver and tolerances
using DifferentialEquations
function r_of_collapse(u, vmin, vmax, M, ρ)
    rvmin(u) = BigFloat(u)/BigFloat(2)     
    m(v) = (BigFloat(M)/BigFloat(2))*(BigFloat(1) + tanh(BigFloat(ρ)*BigFloat(v)))
    B    = - BigFloat(1)/BigFloat(2)
    frhs = (r,p,v) -> -B*(BigFloat(1) - BigFloat(2)*m(v)/r)

    prob = ODEProblem(frhs, rvmin(u), (BigFloat(vmin), BigFloat(vmax)))
    @time begin
    sol  = DifferentialEquations.solve(prob, 
                                       RK4(),
                                       dv=BigFloat(1)/BigFloat(2^12))
    end
    return sol
end

