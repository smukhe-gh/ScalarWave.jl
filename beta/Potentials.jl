#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 04-2018
#--------------------------------------------------------------------

function superposedgaussian{T<:Float64}(u::T, v::T, tarr::Array{T,1}, xarr::Array{T,1}, width::T)::T
    sum = 0.0
    for index in 1:length(xarr)
        u0 = tarr[index] - xarr[index]
        v0 = tarr[index] + xarr[index]
        sum = sum + 10*exp(-(u-u0)^2/width^2)*exp(-(v-v0)^2/width^2)
    end
    return sum
end

function perimeterpotential(Nx::Int, Ny::Int, gsnw::Float64)::Array{Float64,2}
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

    uvgrid = Float64[superposedgaussian(u,v, tarr, xarr, gsnw) for u in chebgrid(6Nx+5), v in chebgrid(6Ny+5)]
    return uvgrid 
end

function piedistribute{T<:Integer}(fbndr::Function, fbndc::Function, frhs::Function, Nx::T, Ny::T, M::T)::Dict{Array{Int,1}, Patch}
    dop = derivOP(Nx, Ny)/M^4
    bop = boundaryOP(Nx, Ny)
    dbase  = Dict{Array{Int, 1}, Patch}()
    for i in 2:2M, k in i-min(i-1,M):min(i-1,M) 
        loc  = [k, i-k]
        v    = -10*perimeterpotential(Int((Nx-5)/6), Int((Ny-5)/6), 0.08)
        rhs  = Float64[chebw(i, Nx)*chebw(i, Ny)*v[i,j] for i in 1:Nx+1, j in 1:Ny+1] 
        bndx = (loc[2]==1) ? (getPatchIC(fbndr, 0, Nx, M, loc[1])) : (getPatchBnd(dbase[loc-[0,1]], 0))
        bndy = (loc[1]==1) ? (getPatchIC(fbndc, 1, Ny, M, loc[2])) : (getPatchBnd(dbase[loc-[1,0]], 1))
        dbase[loc] = calcPatch(bndx, bndy, rhs, dop, bop, loc)
    end
    return dbase
end

