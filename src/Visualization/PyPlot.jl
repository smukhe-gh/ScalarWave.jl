#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2019
# Overload plotting routines from PyPlot
#--------------------------------------------------------------------
using PyPlot
export plot, pcolormesh

function PyPlot. plot(u::Field{S}; plotstyle="-o") where {S}
    x = Field(u.space, x->x)
    plot(x.value,  u.value, plotstyle) 
    return 0
end

function PyPlot. pcolormesh(f::Field{ProductSpace{S1, S2}}) where {S1, S2} 
    u  = Field(f.space.S1, u->u)
    v  = Field(f.space.S2, v->v)
    wu = [integral(f.space.S1, i) for i in 1:size(f.space.S1)]
    wv = [integral(f.space.S2, j) for j in 1:size(f.space.S2)]
    cu = u.value .- (wu/2) 
    cv = v.value .- (wv/2) 
    append!(cu, u.value[end] + (wu[end]/2))
    append!(cv, v.value[end] + (wv[end]/2))
    pcolormesh(cu, cv, f.value, snap=true)
    colorbar()
    return 0
end

