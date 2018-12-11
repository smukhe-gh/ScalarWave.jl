#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 10-2018
# Test the ODE solver with Julia defined types
#--------------------------------------------------------------------

using DifferentialEquations
using DoubleFloats

# using Float64
@time begin
f = (u,p,t) -> (1.01*u)
prob_ode_linear = ODEProblem(f,1/2,(0.0,1.0));
sol = DifferentialEquations.solve(prob_ode_linear, RK4(), dt=1/2^(6))
println(sol)
end

# using DoubleFloats
@time begin
const dfalpha = Double64(101)/Double64(100)
f = (u,p,t) -> (dfalpha*u)
prob_ode_doublefloatlinear = ODEProblem(f, Double64(1)/Double64(2), (Double64(0.0), Double64(1.0)))
sol = DifferentialEquations.solve(prob_ode_doublefloatlinear, RK4(), dt=Double64(1)/Double64(2^6))
println(string(sol[end]))
end

# using BigFloats
@time begin
const bfalpha = BigFloat(101)/BigFloat(100)
f = (u,p,t) -> (bfalpha*u)
prob_ode_bigfloatlinear = ODEProblem(f, BigFloat(1)/BigFloat(2), (BigFloat(0.0), BigFloat(1.0)))
sol = DifferentialEquations.solve(prob_ode_bigfloatlinear, RK4(), dt=BigFloat(1)/BigFloat(2^6))
println(string(sol[end]))
end
