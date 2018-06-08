#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 04-2018
#--------------------------------------------------------------------

function drawPI(Nx::Int, Ny::Int)
    txgrid   = zeros(6*(Nx+1), 6*(Ny+1))
    tcoords  = linspace(-1, 1, 6(Nx+1))  
    xcoords  = linspace(-1, 1, 6(Ny+1))  
    (NX, NY) = (Nx+1, Ny+1) 

    # draw P
    txgrid[2NX:4NX,2NY] = 1
    txgrid[2NX:3NX,3NY] = 1
    
    txgrid[2NX,2NY:3NY] = 1
    txgrid[3NX,2NY:3NY] = 1
    
    # draw I 
    txgrid[2NX:4NX,4NY] = 1
    txgrid[2NX, 4NY - Int(floor(NY/2)):4NY + Int(floor(NY/2))] = 1
    txgrid[4NX, 4NY - Int(floor(NY/2)):4NY + Int(floor(NY/2))] = 1
    
    # flip grid to make PI 'look right' on the UV grid
    txgrid  = flipdim(txgrid, 2)
    
    # store the physical coordinate locations of the points
    tarr = Float64[]
    xarr = Float64[]

    for index in CartesianRange(size(txgrid))
        if txgrid[index] == 1
            push!(tarr, tcoords[index[1]])
            push!(xarr, xcoords[index[2]])
        end
    end
    return (tarr, xarr) 
end

function PIpotential{T<:Float64}(u::T, v::T, width::T, Nx::Int, Ny::Int)::T
    sum = 0.0
    (tarr, xarr) = drawPI(Nx, Ny)
    for index in 1:length(xarr)
        u0 = tarr[index] - xarr[index]
        v0 = tarr[index] + xarr[index]
        sum = sum + exp(-(u-u0)^2/width^2)*exp(-(v-v0)^2/width^2)
    end
    return sum
end

function PIcomputeonPatch(fn::Function, Nx::Int, Ny::Int, M::Int, loc::Array{Int,1})::Array{Float64,2}
    fNpatch = Float64[fn(x,y) for x in chebgrid(Nx, M, loc[1]), y in chebgrid(Ny, M, loc[2])]
    return fNpatch  
end

function PIRHS{T<:Int}(fn::Function, Nx::T, Ny::T, M::T, loc::Array{Int,1})::Array{Float64,2}
    rhs = PIcomputeonPatch(fn, Nx, Ny, M, loc)
    for index in CartesianRange(size(rhs))
        i = index.I[1]
        j = index.I[2]
        rhs[i,j] = chebw(i,Nx)*chebw(j,Ny)*rhs[i,j]
    end
    return rhs
end

function PIdistribute{T<:Integer}(fbndr::Function, fbndc::Function, frhs::Function, Nx::T, Ny::T, M::T)::Dict{Array{Int,1}, Patch}
    dop = derivOP(Nx, Ny)
    bop = boundaryOP(Nx, Ny)
    dbase  = Dict{Array{Int, 1}, Patch}()
    for i in 2:2M, k in i-min(i-1,M):min(i-1,M) 
        loc  = [k, i-k]
        rhs  = PIRHS(frhs, Nx, Ny, M, loc) #Float64[frhs(x,y) for x in chebgrid(Nx, M, loc[1]), y in chebgrid(Ny, M, loc[2])]
        bndx = (loc[2]==1) ? (getPatchIC(fbndr, 0, Nx, M, loc[1])) : (getPatchBnd(dbase[loc-[0,1]], 0))
        bndy = (loc[1]==1) ? (getPatchIC(fbndc, 1, Ny, M, loc[2])) : (getPatchBnd(dbase[loc-[1,0]], 1))
        dbase[loc] = calcPatch(bndx, bndy, rhs, dop, bop, loc) #Patch(loc, rhs)
    end
    return dbase
end
