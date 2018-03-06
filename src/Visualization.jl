#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 03-2018
#--------------------------------------------------------------------

function setcolormap(vec::Array{Float64,1}, map::String, samples::Int)
    nvec = (vec .- minimum(vec))./(maximum(vec) .- minimum(vec))
    clrs = colormap(map, samples+1)
    return clrs[round.(Int, (nvec*samples)+1)]
end

function drawpatch(patch::Array{Float64,2})
    (Nx, Ny) = size(patch) .- 1
    (wx, wy) = (chebweights(Nx), chebweights(Ny))
    (lx, ly) = (sort(2.0 - cumsum(wx))*200, sort(2.0 - cumsum(wy))*200)
    cmap     = reshape(setcolormap(vec(patch), "Blues", 100000), size(patch))
    canvas   = Drawing(800, 800, "luxor-patch.pdf")

    #-----------------------------------------------
    # draw patch
    #-----------------------------------------------
    origin(400, 650)
    rotate(-3pi/4)
    setline(0.4)
    for i in 2:Nx, j in 2:Ny
       sethue((cmap[i,j].r, cmap[i,j].g, cmap[i,j].b))
       rect(lx[i], ly[j], lx[i+1] - lx[i], ly[j+1] - ly[j], :fill)
    end

    #-----------------------------------------------
    # set ticks
    #-----------------------------------------------
    sethue("black")
    setline(0.4)
    (lx0, ly0) = (lx[2], ly[2])
    (lxe, lye) = (lx[end], ly[end])
    (lxm, lym) = (lx0, ly0)./2 .+ (lxe, lye)./2

    # corners
    line(Point(lx0, ly0), Point(lx0, ly0-10), :stroke)
    line(Point(lx0, ly0), Point(lx0-10, ly0), :stroke)
    line(Point(lxe, lye), Point(lxe, lye+10), :stroke)
    line(Point(lxe, lye), Point(lxe+10, lye), :stroke)
    line(Point(lxe, ly0), Point(lxe, ly0-10), :stroke)
    line(Point(lxe, ly0), Point(lxe+10, ly0), :stroke)
    line(Point(lx0, lye), Point(lx0-10, lye), :stroke)
    line(Point(lx0, lye), Point(lx0, lye+10), :stroke)
    line(Point(lx0, ly0), Point(lx0, ly0-10), :stroke)
    line(Point(lx0, ly0), Point(lx0-10, ly0), :stroke)

    # mid points
    line(Point(lxm, ly0), Point(lxm, ly0-10), :stroke)
    line(Point(lx0, lym), Point(lx0-10, lym), :stroke)
    line(Point(lxm, lye), Point(lxm, lye+10), :stroke)
    line(Point(lxe, lym), Point(lxe+10, lym), :stroke)

    # edges
    line(Point(lx0, ly0), Point(lx0, lye), :stroke)
    line(Point(lx0, ly0), Point(lxe, ly0), :stroke)
    line(Point(lxe, ly0), Point(lxe, lye), :stroke)
    line(Point(lx0, lye), Point(lxe, lye), :stroke)

    #-----------------------------------------------
    # tick labels
    #-----------------------------------------------
    settext("(1, 1)", Point(lx0 - 18, ly0 - 18);
            halign = "center",
            valign = "top")
    settext("(-1, -1)", Point(lxe + 18, lye + 18);
            halign = "center",
            valign = "bottom")
    settext("(-1, 1)", Point(lx0 - 20, lye + 20);
            halign = "center",
            valign = "center")
    settext("(1, -1)", Point(lxe + 20, ly0 - 20);
            halign = "center",
            valign = "center")

    #-----------------------------------------------
    # colorbar
    #-----------------------------------------------
    rotate(3pi/4)
    x   = collect(linspace(-200, 200, 100))
    clr = colormap("Blues", 100)

    for i in 1:99
        sethue((clr[i].r, clr[i].g, clr[i].b))
        rect(x[i], 92.5, x[i+1] - x[i], 15, :fill)
    end

    #-----------------------------------------------
    # colobar ticks
    #-----------------------------------------------
    sethue("black")
    setline(0.4)
    line(Point(-200, 107.5), Point(-200, 85), :stroke)
    line(Point(0, 107.5), Point(0, 85), :stroke)
    line(Point(200, 107.5), Point(200, 85), :stroke)
    min = round(Int, findmin(patch)[1])
    max = round(Int, findmax(patch)[1])
    med = round(Int, (min + max)/2)

    settext("$min", Point(-200, 80);
                halign = "center",
                valign = "center")
    settext("$med", Point(0, 80);
                halign = "center",
                valign = "center")
    settext("$max", Point(200, 80);
                halign = "center",
                valign = "center")

    #-----------------------------------------------
    #arrow
    #-----------------------------------------------
    sethue("black")
    setline(0.4)
    origin(700, 100)
    arrow(O, O .+ Point(40, -40))
    settext("u", O .+ Point(40, -40);
                halign = "top",
                valign = "right")
    rotate(-pi/2)
    arrow(O, O .+ Point(40, -40))
    settext("v", O .+ Point(40, -40);
                halign = "bottom",
                valign = "right")

    finish()
    return canvas
end
