#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------
# Functions to test a generic future implementation.
# Doen not test the functions in src 

@test typeof(@spawn 2) === Future

function pascal!(A::Array{Float64,2})::Array{Float64,2}
    N = size(A)[1]
    for i in 2:N, j in 2:N
        A[i,j] = A[i,j-1] + A[i-1,j]
    end
    return A
end

function fpascal(A::Future)::Future
    return @spawn pascal!(fetch(A))
end

