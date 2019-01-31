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

# plot as an image
# imshow(ϕcart.value, origin="lower", extent=[-5, 15, -5, 5])

# now plot a contour
CV = Taylor{V, 200, 15//1, -5//1}
CU = Taylor{U, 100,  5//1, -5//1}
UU = Field(CU, x->x)
VV = Field(CV, x->x)

p1 = contour(VV.value, UU.value, ϕcart.value)
xlabel(L"V")
ylabel(L"U")
#show()
close()


#--------------------------------------------------------------------
# Now make 2D plots first
#--------------------------------------------------------------------

struct R end
struct T end

M    =  1.0
Rmax =  3.0
Rmin =  1.0
Tmax =  1.0
Tmin = -1.0

SR = GaussLobatto{R, 20, Rmax, Rmin}
ST = GaussLobatto{T, 20, Tmax, Tmin}
S  = ProductSpace{SR, ST}

u = Field(S, (T,R)->find_U_of_tr(T,R,M))
v = Field(S, (T,R)->find_V_of_tr(T,R,M))

# Now plot
CSR = Taylor{R, 12, rationalize(Rmax),  rationalize(Rmin)}
CST = Taylor{T, 10, rationalize(Tmax),  rationalize(Tmin)}
CS  = ProductSpace{CSR, CST}

gt = Field(CST, t->t)
gr = Field(CSR, r->r)
uu = Field(CS,  (t,r)->u(t,r))
vv = Field(CS,  (t,r)->v(t,r))
tt = Field(CS,  (t,r)->t)
rr = Field(CS,  (t,r)->r)

outgoing = contour(gr.value, gt.value, uu.value, levels=12)
# incoming = contour(gr.value, gt.value, vv.value, levels=20)
xlabel(L"R")
ylabel(L"T")
colorbar()
show()
