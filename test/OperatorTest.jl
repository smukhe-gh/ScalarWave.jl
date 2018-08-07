
#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
# Test summation-by-parts property for our derivative operator
# See Summation by Parts, Projections, and Stability. I Olsson 1995
#--------------------------------------------------------------------

D  = [-1.0  1.0  0.0  0.0;
      -0.5  0.0  0.5  0.0;
       0.0 -0.5  0.0  0.5;
       0.0  0.0 -1.0  1.0]

W  = [ 0.5  0.0  0.0  0.0;
       0.0  1.0  0.0  0.0;
       0.0  0.0  1.0  0.0;
       0.0  0.0  0.0  0.5]

specD = Float64[chebd(i,j, 3) for i in 1:4, j in 1:4]
specW = diagm(Float64[chebw(i, 3) for i in 1:4])

u = [x^2 for x in linspace(-1, 1, 4)]
v = [x   for x in linspace(-1, 1, 4)]

@test u'*W*D*v + v'*W*D*u == u[end]*v[end] - u[1]*v[1]
@test_broken u'*specW*specD*v + v'*specW*specD*u == u[end]*v[end] - u[1]*v[1]

