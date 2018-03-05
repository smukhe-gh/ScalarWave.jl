#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 03-2018
#--------------------------------------------------------------------

function generatecmap(patch::Array{Float64,2}, N::Int)
    npatch = (patch .- minimum(patch))./(maximum(patch) .- minimum(patch))
    colors = colormap("Blues", N+1)
    return reshape(colors[vec(round.(Int, npatch*N+1))], size(npatch))
end

function drawpatch(patch::Array{Float64, 2})

    (Nx, Ny) = size(patch) .- 1
    (wx, wy) = (chebweights(Nx), chebweights(Ny))
    (lx, ly) = (sort(2.0 - cumsum(wx))*200, sort(2.0 - cumsum(wy))*200)
    cmap = generatecmap(patch, 100000)
    canvas = Drawing(800, 800, "plotpatch.pdf")

    # draw patch
    origin(400, 650)
    rotate(-3pi/4)
    setline(0.4)
    for i in 2:Nx, j in 2:Ny
       sethue((cmap[i,j].r, cmap[i,j].g, cmap[i,j].b))
       rect(lx[i], ly[j], lx[i+1] - lx[i], ly[j+1] - ly[j], :fill)
    end

    # set ticks
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

    # tick labels
    settext("(1, 1)", Point(lx0 - 18, ly0 - 18);
            halign = "center",
            valign = "top",
            markup = false)

    settext("(-1, -1)", Point(lxe + 18, lye + 18);
            halign = "center",
            valign = "bottom",
            markup = false)
    settext("(-1, 1)", Point(lx0 - 20, lye + 20);
            halign = "center",
            valign = "center",
            markup = false)
    settext("(1, -1)", Point(lxe + 20, ly0 - 20);
            halign = "center",
            valign = "center",
            markup = false)

    # colorbar
    rotate(3pi/4)
    # line(Point(-200, +100), Point(200, +100))
    # strokepath()

    x   = collect(linspace(-200, 200, 100))
    clr = colormap("Blues", 100)

    for i in 1:99
        sethue((clr[i].r, clr[i].g, clr[i].b))
        rect(x[i], 92.5, x[i+1] - x[i], 15, :fill)
    end

    # colobar ticks
    sethue("black")
    setline(0.4)
    line(Point(-200, 107.5), Point(-200, 85), :stroke)
    line(Point(200, 107.5), Point(200, 85), :stroke)
    line(Point(0, 107.5), Point(0, 85), :stroke)

    min = round(Int, findmin(patch)[1])
    max = round(Int, findmax(patch)[1])
    med = round(Int, (min + max)/2)

    settext("$min", Point(-200, 80);
                halign = "center",
                valign = "center",
                markup = false)
    settext("$med", Point(0, 80);
                halign = "center",
                valign = "center",
                markup = false)
    settext("$max", Point(200, 80);
                halign = "center",
                valign = "center",
                markup = false)

    #arrow
    origin(700, 100)
    arrow(O, O .+ Point(40, -40))
    settext("u", O .+ Point(40, -40);
                halign = "top",
                valign = "right",
                markup = false)
    rotate(-pi/2)
    arrow(O, O .+ Point(40, -40))
    settext("v", O .+ Point(40, -40);
                halign = "bottom",
                valign = "right",
                angle  = 0,
                markup = false)
    finish()
    return canvas
end
