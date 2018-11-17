using Plots, LinearAlgebra, ApproxFun

dx = -1.0 .. 1.0; dt = 2.0 .. 3.0
d = dx × dt
Dx = ApproxFun.Derivative(d,[1,0]); Dt = ApproxFun.Derivative(d,[0,1])
u0 = Fun(x->exp(-20x^2),dx)
u0 = Fun((x,_) -> exp(-20x^2), ApproxFun.Chebyshev(-1.0..1.0)⊗ConstantSpace(ApproxFun.Point(2.0)))

@time u = \([I⊗ldirichlet(dt); ldirichlet(dx)⊗I; Dt+Dx], [u0; 0; 0];
            tolerance=1E-3)
pyplot()
plot(u)
savefig("../output/convection-u.pdf")
