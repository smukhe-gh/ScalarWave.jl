#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

@btime fdistribute(x->sin(x), y->sin(y), (x,y)->0, 40, 40, 20)
@btime distribute(x->sin(x), y->sin(y), (x,y)->0,  40, 40, 20)
