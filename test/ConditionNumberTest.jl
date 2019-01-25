#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 04-2018
# Check how condition numbers grow with resolution 
#--------------------------------------------------------------------

using LinearAlgebra, Plots

N = 10
condition_number = zeros(N)
for n in 1:N
    SUV = ProductSpace{GaussLobatto(V,n),
                       GaussLobatto(U,n)}
    𝔻𝕍, 𝔻𝕌 = derivative(SUV) 
    𝔹 = boundary(Null, SUV) 
    condition_number[n] = cond((1/n^2)*(𝔻𝕌*𝔻𝕍) + 𝔹)
end

display(condition_number)
println("")
