#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 04-2018
#--------------------------------------------------------------------

patch = CompactifiedMinkowskidistribute(u->sin(pi*u), v->sin(pi*v), (u,v)->0, 20, 20)
drawmultipatch(patch)
