#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2019
# Overload plotting routines from PyPlot
#--------------------------------------------------------------------

function PyPlot. plot(u::Field{S}; collocation=false,
                                   plotstyle="r-", 
                                   npoints=100) where {S<:GaussLobatto{Tag, N, max, min}} where {Tag, N, max, min}
    
    cc = basistransform(u)
    SS = Taylor{Tag, npoints, max, min}
    uu = Field(SS, x->cc(x))
    xx = Field(SS, x->x)

    if collocation == false
        return plot(xx.value, uu.value, plotstyle)
    else
        x = Field(S, x->x)
        plotstyle="o"
        return plot(x.value,  u.value, plotstyle) 
    end
end

function PyPlot. pcolormesh(u::Field{S}) where {S<:ProductSpace{GaussLobatto{Tag1, N1, max1, min1},
                                                                GaussLobatto{Tag2, N2, max2, min2}}} where {Tag1, N1, max1, min1,
                                                                                                            Tag2, N2, max2, min2}

    # compute the coordinates of the corner
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

function PyPlot. contourf(u::Field{S}, fill=true) where {S<:ProductSpace{GaussLobatto{Tag1, N1, max1, min1},
                                                                         GaussLobatto{Tag2, N2, max2, min2}}} where {Tag1, N1, max1, min1,
                                                                                                                     Tag2, N2, max2, min2}
    cc  = basistransform(u)
    xx1 = Field(Taylor{Tag1, N1, max1, min1}, x->x)
    xx2 = Field(Taylor{Tag2, N2, max2, min2}, x->x)
    uu  = Field(ProductSpace{Taylor{Tag1, N1, max1, min1},
                             Taylor{Tag2, N2, max2, min2}}, (x,y)->cc(x,y)) 

    if fill==true
        return contourf(xx1.value, xx2.value, uu.value)
    else
        return contour(xx1.value, xx2.value, uu.value)
    end
end


