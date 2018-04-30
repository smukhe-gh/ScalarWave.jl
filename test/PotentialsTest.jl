#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 04-2018
#--------------------------------------------------------------------

patch = piedistribute(x->0, y->10*sin(pi*y), (x,y)->0, 20, 20, 1)
drawpatch(patch[[1,1]].value)

