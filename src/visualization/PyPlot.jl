#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2019
# Overload plotting routines from PyPlot
#--------------------------------------------------------------------
# FIXME: Fix type inference to prevent type conversion in Taylor basis.

function PyPlot. plot(u::Field{S}; plotstyle="o") where {S<:GaussLobatto{Tag, N, max, min}} where {Tag, N, max, min}
    x = Field(S, x->x)
    plot(x.value,  u.value, plotstyle) 
    return 0
end

function PyPlot. plot(u::Field{S}, npoints::Int; plotstyle="-") where {S<:GaussLobatto{Tag, N, max, min}} where {Tag, N, max, min}
    cc = basistransform(u)
    SS = Taylor{Tag, npoints, max, min}
    uu = Field(SS, x->cc(x))
    xx = Field(SS, x->x)
    plot(xx.value, uu.value, plotstyle)
    return 0
end

function PyPlot. pcolormesh(u::Field{S}) where {S<:ProductSpace{GaussLobatto{Tag1, N1, max1, min1},
                                                                GaussLobatto{Tag2, N2, max2, min2}}} where {Tag1, N1, max1, min1,
                                                                                                            Tag2, N2, max2, min2}
    xx1 = Field(GaussLobatto{Tag1, N1, max1, min1}, x->x)
    xx2 = Field(GaussLobatto{Tag2, N2, max2, min2}, x->x)
    w1  = [chebw(i, N1) for i in 1:N1+1]*((max1 - min1)/2)
    w2  = [chebw(i, N2) for i in 1:N2+1]*((max2 - min2)/2)
    c1  = xx1.value .- (w1/2) 
    c2  = xx2.value .- (w2/2) 
    append!(c1, xx1.value[end] + (w1[end]/2))
    append!(c2, xx2.value[end] + (w2[end]/2))
    pcolormesh(c1, c2, u.value, snap=true)
    return 0
end

function PyPlot. pcolormesh(u::Field{S}, npoints::Int; globalmax=nothing, globalmin=nothing) where {S<:ProductSpace{GaussLobatto{Tag1, N1, max1, min1},
                                                                             GaussLobatto{Tag2, N2, max2, min2}}} where {Tag1, N1, max1, min1,
                                                                                                                         Tag2, N2, max2, min2}
    xx1 = Field(Taylor{Tag1, npoints, max1, min1}, x->x)
    xx2 = Field(Taylor{Tag2, npoints, max2, min2}, x->x)
    w1  = [1 for i in 1:npoints+1]*((max1 - min1)/2)
    w2  = [1 for i in 1:npoints+1]*((max2 - min2)/2)
    c1  = xx1.value .- (w1/2) 
    c2  = xx2.value .- (w2/2) 
    append!(c1, xx1.value[end] + (w1[end]/2))
    append!(c2, xx2.value[end] + (w2[end]/2))

    cc = basistransform(u)
    uu = Field(ProductSpace{Taylor{Tag1, npoints, max1, min1},
                            Taylor{Tag2, npoints, max2, min2}}, (x,y)->cc(x,y))

    # N.B: Unnecessary type conversion. Need to fix this.
    pcolormesh(Float64.(c1), Float64.(c2), Float64.(uu.value), snap=true, vmin=globalmin, vmax=globalmax)        
    return 0
end

function PyPlot. contourf(u::Field{S}, npoints::Int; globalmax=nothing, globalmin=nothing, globallevels=nothing) where {S<:ProductSpace{GaussLobatto{Tag1, N1, max1, min1},
                                                                           GaussLobatto{Tag2, N2, max2, min2}}} where {Tag1, N1, max1, min1,
                                                                                                                       Tag2, N2, max2, min2}
   cc  = basistransform(u)
   xx1 = Field(Taylor{Tag1, npoints, max1, min1}, x->x)
   xx2 = Field(Taylor{Tag2, npoints, max2, min2}, x->x)
   uu  = Field(ProductSpace{Taylor{Tag1, npoints, max1, min1},
                            Taylor{Tag2, npoints, max2, min2}}, (x,y)->cc(x,y)) 
   contourf(Float64.(xx1.value), Float64.(xx2.value), Float64.(uu.value), vmin=globalmin, vmax=globalmax, levels=globallevels)
   return 0
end


function PyPlot. contour(u::Field{S}, npoints::Int; globalmax=nothing, globalmin=nothing, globallevels=nothing) where {S<:ProductSpace{GaussLobatto{Tag1, N1, max1, min1},
                                                                           GaussLobatto{Tag2, N2, max2, min2}}} where {Tag1, N1, max1, min1,
                                                                                                                       Tag2, N2, max2, min2}
   cc  = basistransform(u)
   xx1 = Field(Taylor{Tag1, npoints, max1, min1}, x->x)
   xx2 = Field(Taylor{Tag2, npoints, max2, min2}, x->x)
   uu  = Field(ProductSpace{Taylor{Tag1, npoints, max1, min1},
                            Taylor{Tag2, npoints, max2, min2}}, (x,y)->cc(x,y)) 
   contour(Float64.(xx1.value), Float64.(xx2.value), Float64.(uu.value), vmin=globalmin, vmax=globalmax, levels=globallevels)
   return 0
end

function PyPlot. plot(nest::Dict, npoints::Int; plotstyle="-")
    for (loc, value) in nest 
        plot(value, npoints)
    end
    return 0
end

function PyPlot. contourf(nest::Dict{Array{Int64,1}, Union{Field, Dict}}, npoints::Int64; globalmax=nothing, globalmin=nothing, globallevels=nothing)
    for (loc, value) in nest 
        contourf(value, npoints, globalmax=globalmax, globalmin=globalmin, globallevels=globallevels)
    end
    return 0
end

function PyPlot. contour(nest::Dict{Array{Int64,1}, Union{Field, Dict}}, npoints::Int64; globalmax=nothing, globalmin=nothing, globallevels=nothing)
    for (loc, value) in nest
        contour(value, npoints, globalmax=globalmax, globalmin=globalmin, globallevels=globallevels)
    end
    return 0
end

function Base. maximum(nest::Dict{Array{Int64,1}, Union{Field, Dict}})
    localmax = []
    for (loc, value) in nest
        append!(localmax, maximum(value))
    end
    return maximum(localmax)
end

function Base. minimum(nest::Dict{Array{Int64,1}, Union{Field, Dict}})
    localmin = []
    for (loc, value) in nest
        append!(localmin, minimum(value))
    end
    return minimum(localmin)
end

function levels(nest::Dict{Array{Int64,1}, Union{Field, Dict}}; globallength=20)
    return collect(range(minimum(nest), stop=maximum(nest), length=globallength))
end
