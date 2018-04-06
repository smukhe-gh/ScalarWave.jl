#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

fbase = fdistribute(x->sin(x), y->sin(y), (x,y)->0, 2, 2, 2)
dbase = distribute(x->sin(x), y->sin(y), (x,y)->0, 2, 2, 2)

@test fbase[[1,1]].value ≈ dbase[[1,1]].value
@test fbase[[1,2]].value ≈ dbase[[1,2]].value
@test fbase[[2,1]].value ≈ dbase[[2,1]].value
@test fbase[[2,2]].value ≈ dbase[[2,2]].value
