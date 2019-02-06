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
# plot(coords(yL).value, yL.value, "r")
# plot(coords(yR).value, yR.value, "b")

# now join the space 
yLR = coarsen(yL, yR)
# plot(coords(yLR).value, yLR.value, "go")

# show()

#--------------------------------------------------------------------
# Now do a 2D refinement 
#--------------------------------------------------------------------

SUV = ProductSpace{GaussLobatto{V, 20,  5.0, -5.0},
                   GaussLobatto{U, 20,  5.0, -5.0}}

ϕ   = Field(SUV, (U,V)->sinpi(U)*cospi(V))

SUV = ProductSpace{GaussLobatto{V, 20,  15.0, 5.0},
                   GaussLobatto{U, 20,  15.0, 5.0}}

ψ   = Field(SUV, (U,V)->sinpi(U)*cospi(V))

function PyPlot. imshow(u::Field{S}) where {S<:ProductSpace{GaussLobatto{TagV, NV, maxV, minV}, 
                                                            GaussLobatto{TagU, NU, maxU, minU}}} where {TagV, NV, maxV, minV, 
                                                                                                        TagU, NU, maxU, minU}
    ucoefficents = basistransform(u)
    Scartesian   = ProductSpace{Taylor{TagV, 100, rationalize(maxV), rationalize(minV)},
                                Taylor{TagU, 100, rationalize(maxU), rationalize(minU)}}
    ucartesian   = Field(Scartesian, (x,y)->ucoefficents(x,y))
    return imshow(ucartesian.value, origin="lower", extent=[minV, maxV, minU, maxU])
end

p = imshow(ϕ)
q = imshow(ψ)
show()
