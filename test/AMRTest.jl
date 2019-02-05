#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2019
# Plot fields using PyPlot
#--------------------------------------------------------------------

using PyPlot

#--------------------------------------------------------------------
# Refine in 1D first
#--------------------------------------------------------------------

# construct the field
S = GaussLobatto{U, 100, 3.0, -3.0}
x = Field(S, x->x)
y = sin(x)

# divide the space in 2
yL, yR = refine(y)

function coords(u::Field{S}) where {S<:GaussLobatto{Tag, N, max, min}} where {Tag, N, max, min}
    coords = zeros(N+1)
    for i in 1:N+1
        coords[i] = collocation(S, i)
    end
    return Field(S, coords)
end

# pass it to the plotting routine
plot(coords(yL).value, yL.value, "r")
plot(coords(yR).value, yR.value, "b")

# now join the space 
yLR = coarsen(yL, yR)
plot(coords(yLR).value, yLR.value, "go")

show()
