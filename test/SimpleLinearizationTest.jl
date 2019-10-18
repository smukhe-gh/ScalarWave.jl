#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 10-2019
# Test linearization
#--------------------------------------------------------------------

function f!(F, x)
    F[1] = (x[1]+3)*(x[2]^3-7)+18
    F[2] = sin(x[2]*exp(x[1])-1)
end


function j!(J, x)
    J[1, 1] = x[2]^3-7
    J[1, 2] = 3*x[2]^2*(x[1]+3)
    u = exp(x[1])*cos(x[2]*exp(x[1])-1)
    J[2, 1] = x[2]*u
    J[2, 2] = u
end

F = zeros(2)
J = zeros(2,2)
x = rand(2)

f!(F,x)
j!(J,x)
@show F
@show J

using ForwardDiff
@show ForwardDiff.jacobian(f!, F, x)

