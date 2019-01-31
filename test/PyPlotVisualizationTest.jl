#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2019
# Plot fields using PyPlot
#--------------------------------------------------------------------

using PyPlot

#--------------------------------------------------------------------
# Make 1D plots first
#--------------------------------------------------------------------

# construct the field
S = GaussLobatto{U, 10, 3.0, -3.0}
x = Field(S, x->x)
y = sin(x)

# now interpolate
xinterp = range(minimum(S), stop=maximum(S), length=1000)
yinterp = y.(xinterp)

# Can we do this differently?
SC = Taylor{U, 1000, 3//1, -3//1}
xtaylorinterp = Field(SC, x->x)
ytaylorinterp = Field(SC, x->y(x)) 

# pass it to the plotting routine
plot(xtaylorinterp.value, ytaylorinterp.value, "b-")
plot(xinterp, yinterp)
plot(x.value, y.value, "ro")
# show()
close()
#--------------------------------------------------------------------
# Now make 2D plots first
#--------------------------------------------------------------------

SUV = ProductSpace{GaussLobatto{V, 10, 15//1, -5//1},
                   GaussLobatto{U, 14,  5//1, -5//1}}

ϕ   = Field(SUV, (U,V)->exp(-U^2) + exp(-V^2))

CUV = ProductSpace{Taylor{V, 200, 15//1, -5//1},
                   Taylor{U, 100,  5//1, -5//1}}


ϕcart = Field(CUV, (U,V)->exp(-U^2) + exp(-V^2))
imshow(ϕcart.value, origin="lower", extent=[-5, 15, -5, 5])
xlabel(L"V")
ylabel(L"U")
show()
