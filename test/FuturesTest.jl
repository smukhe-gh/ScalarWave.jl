#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

fdbase = fdistribute(x->sin(x), y->sin(y), (x,y)->0, 2, 2, 2)
dbase  = distribute(x->sin(x), y->sin(y), (x,y)->0, 2, 2, 2)
@test fetch(dbase[[1,1]]) == dbase[[1,1]]
@test fetch(dbase[[1,2]]) == dbase[[1,2]]
@test fetch(dbase[[2,1]]) == dbase[[2,1]]
@test fetch(dbase[[2,2]]) == dbase[[2,2]]
