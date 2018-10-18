#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 03-2018
#--------------------------------------------------------------------

function drawtensorfield(g::Metric, filename)
    eigU, eigV = eigvals(g)

    # XXX: Hack to minimize changes in the rest of the code
    dbase = Dict([1,1]=>Patch([1,1], eigU.value))

    M   = 1
    AP  = Float64[]

    #-----------------------------------------------
    # Create color pallete
    #-----------------------------------------------
    for mx in 1:M, my in 1:M
        AP = vcat(AP, vec(dbase[[mx, my]].value))
    end
    globalcmap = setcolormap(AP, "Blues", length(AP))
   

    #-----------------------------------------------
    # Patch sizing 
    #-----------------------------------------------
    samplepatch = dbase[[1,1]].value
    (Nx, Ny) = size(samplepatch) .- 1
    (wx, wy) = (chebweights(Nx), chebweights(Ny))
    (lx, ly) = (sort(2.0 - cumsum(wx))*200, sort(2.0 - cumsum(wy))*200)

    #-----------------------------------------------
    # set-up canvas
    #-----------------------------------------------
    #canvas   = Drawing(800, 800, "$filename.pdf")
    canvas   = Drawing(800, 800, "$filename.pdf")
    push!(lx, 400.0)
    push!(ly, 400.0)
    origin(400, 650)
    rotate(-3pi/4)
   
    setline(0.4)
    for i in 1:Nx+1, j in 1:Ny+1
        rect(lx[i], ly[j], lx[i+1] - lx[i], ly[j+1] - ly[j], :stroke)
    end

    #-----------------------------------------------
    # set ticks
    #-----------------------------------------------
    sethue("black")
    setline(0.4)
    (lx0, ly0) = (lx[1], ly[1])
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
    
    gsave()
    rotate(3pi/4)
    x   = collect(linspace(-200, 200, length(AP)))
    clr = colormap("Blues", length(AP))

    for i in 1:length(AP)-1
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
    
    min = findmin(AP)[1] 
    max = findmax(AP)[1]
    med = (min + max)/2

    settext("$min", Point(-200, 80);
                halign = "center",
                valign = "center")
    settext("$med", Point(0, 80);
                halign = "center",
                valign = "center")
    settext("$max", Point(200, 80);
                halign = "center",
                valign = "center")

    grestore()

    
    #-----------------------------------------------
    # draw patches
    #-----------------------------------------------
    
    for mx in 1:M, my in 1:M
        gsave()
        Luxor.translate(((my-1)/2M)*(maximum(lx)+400), ((mx-1)/2M)*(maximum(lx)+400))
        scale(1/M)
        patch     = dbase[[mx, my]].value
        localcmap = globalcmap[(sub2ind((M, M), mx, my)-1)*length(patch) + 1:(sub2ind((M, M), mx, my))*length(patch)] 
        cmap = reshape(localcmap, size(patch))

        for i in 1:Nx+1, j in 1:Ny+1
           sethue((cmap[i,j].r, cmap[i,j].g, cmap[i,j].b))
           ellipse(Point((lx[i] + lx[i+1])/2, (ly[j] + ly[j+1])/2), eigU.value[i,j]*50, eigV.value[i,j]*50, :fill)
           sethue("black")
           ellipse(Point((lx[i] + lx[i+1])/2, (ly[j] + ly[j+1])/2), eigU.value[i,j]*50, eigV.value[i,j]*50, :stroke)
        end
        
        grestore()
        
    end 
    
    #-----------------------------------------------
    #arrow
    #-----------------------------------------------
    sethue("black")
    setline(0.4)
    origin(700, 100)
    Luxor.arrow(O, O .+ Point(40, -40))
    settext("v", O .+ Point(40, -40);
                halign = "top",
                valign = "right")
    rotate(-pi/2)
    Luxor.arrow(O, O .+ Point(40, -40))
    settext("u", O .+ Point(40, -40);
                halign = "bottom",
                valign = "right")
    
    finish()
    return canvas
end
