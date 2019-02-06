#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2019
# Overload plotting routines from PyPlot
#--------------------------------------------------------------------
# FIXME: Fix type inference to prevent type conversion in Taylor basis.

function PyPlot. plot(u::Field{S}; plotstyle="ro") where {S<:GaussLobatto{Tag, N, max, min}} where {Tag, N, max, min}
    x = Field(S, x->x)
    return plot(x.value,  u.value, plotstyle) 
end

function PyPlot. plot(u::Field{S}, npoints::Int; plotstyle="r-") where {S<:GaussLobatto{Tag, N, max, min}} where {Tag, N, max, min}
    cc = basistransform(u)
    SS = Taylor{Tag, npoints, max, min}
    uu = Field(SS, x->cc(x))
    xx = Field(SS, x->x)
    return plot(xx.value, uu.value, plotstyle)
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
    return pcolormesh(c1, c2, u.value, snap=true)
end

function PyPlot. pcolormesh(u::Field{S}, npoints::Int) where {S<:ProductSpace{GaussLobatto{Tag1, N1, max1, min1},
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
    return pcolormesh(Float64.(c1), Float64.(c2), Float64.(uu.value), snap=true)        
end

function PyPlot. contourf(u::Field{S}, npoints::Int) where {S<:ProductSpace{GaussLobatto{Tag1, N1, max1, min1},
                                                                           GaussLobatto{Tag2, N2, max2, min2}}} where {Tag1, N1, max1, min1,
                                                                                                                       Tag2, N2, max2, min2}
   cc  = basistransform(u)
   xx1 = Field(Taylor{Tag1, npoints, max1, min1}, x->x)
   xx2 = Field(Taylor{Tag2, npoints, max2, min2}, x->x)
   uu  = Field(ProductSpace{Taylor{Tag1, npoints, max1, min1},
                            Taylor{Tag2, npoints, max2, min2}}, (x,y)->cc(x,y)) 
   return contour(Float64.(xx1.value), Float64.(xx2.value), Float64.(uu.value))
end


function PyPlot. contour(u::Field{S}, npoints::Int) where {S<:ProductSpace{GaussLobatto{Tag1, N1, max1, min1},
                                                                          GaussLobatto{Tag2, N2, max2, min2}}} where {Tag1, N1, max1, min1,
                                                                                                                      Tag2, N2, max2, min2}
    cc  = basistransform(u)
    xx1 = Field(Taylor{Tag1, npoints, max1, min1}, x->x)
    xx2 = Field(Taylor{Tag2, npoints, max2, min2}, x->x)
    uu  = Field(ProductSpace{Taylor{Tag1, npoints, max1, min1},
                             Taylor{Tag2, npoints, max2, min2}}, (x,y)->cc(x,y)) 
    return contour(Float64.(xx1.value), Float64.(xx2.value), Float64.(uu.value))
end


