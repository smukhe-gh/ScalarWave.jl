#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2019
# Overload plotting routines from PyPlot
#--------------------------------------------------------------------
using PyPlot
export plot, pcolormesh, contour, contourf

function PyPlot. plot(u::Field{S}; plotstyle="-o", label="") where {S}
    x = Field(u.space, x->x)
    plot(x.value,  u.value, plotstyle, label=label) 
end

# FIXME: This plotting routine doesn't work quite right. 
function PyPlot. pcolormesh(f::Field{ProductSpace{S1, S2}}) where {S1, S2} 
    u  = Field(f.space.S1, u->u)
    v  = Field(f.space.S2, v->v)
    wu = [integral(f.space.S1, i) for i in 1:size(f.space.S1)]
    wv = [integral(f.space.S2, j) for j in 1:size(f.space.S2)]
    cu = u.value .- (wu/2) 
    cv = v.value .- (wv/2) 
    append!(cu, u.value[end] + (wu[end]/2))
    append!(cv, v.value[end] + (wv[end]/2))
    pcolormesh(cv, cu, f.value, snap=true)
    xlabel("v")
    ylabel("u")
    colorbar()
end

function PyPlot. contour(f::Field{ProductSpace{S1, S2}}, levels::Union{Array{Number, 1}}, Int) where {S1, S2} 
    u  = Field(f.space.S1, u->u)
    v  = Field(f.space.S2, v->v)
    cp = contour(v.value, u.value, f.value, levels, colors="k")
    clabel(cp, inline=1, fontsize=5, colors="k")
    xlabel("v")
    ylabel("u")
    colorbar()
end

function PyPlot. contourf(f::Field{ProductSpace{S1, S2}}, levels::Union{Array{Number, 1}, Int}) where {S1, S2} 
    u  = Field(f.space.S1, u->u)
    v  = Field(f.space.S2, v->v)
    contourf(v.value, u.value, f.value, levels)
    xlabel("v")
    ylabel("u")
    colorbar()
end

function PyPlot. contourf(grid::Array{Field, 2}, levels::Int)
    gridmin = minimum((minimum.(grid)))
    gridmax = maximum((maximum.(grid)))
    gridlvl = collect(range(gridmin, stop=gridmax, length=levels))
    for index in CartesianIndices(grid) 
        f  = grid[index]
        u  = Field(f.space.S1, u->u)
        v  = Field(f.space.S2, v->v)
        contourf(v.value, u.value, f.value, levels, vmin = gridmin, vmax = gridmax)
    end
    xlabel("v")
    ylabel("u")
    colorbar()
end

