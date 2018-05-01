#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 04-2018
#--------------------------------------------------------------------

patch = PIdistribute(u->0, v->0, (u,v)->-sin(10*pi*u)*sin(10*pi*v)*PIpotential(u,v,0.05, 20, 20), 20, 20, 8)
#patch = PIdistribute(u->0, v->0, (u,v)->-sin(10*pi*(u+v))*PIpotential(u,v,0.05, 40, 40), 40, 40, 8)
drawmultipatch(patch)
