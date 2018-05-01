#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 04-2018
# Scalar Wave equation in compactified Minkowski background in 
# double null coordinates
# See H.12 in Spacetime and Geometry (Carroll)
#--------------------------------------------------------------------

function ginv(U::Float64, V::Float64)::Float64
    return -2*(cos(U)*cos(V))^2
end

function detg(U::Float64, V::Float64)::Float64
    return -1/(2*cos(U)*cos(V))^2
end

function coordtransUV2tr(U::Float64, V::Float64)::(Float64, Float64)
    u = tan(U)
    v = tan(V)
    t = (v + u)/2
    r = (v - u)/2
    return (t,r)
end

function coordtranstr2UV(t::Float64, r::Float64)::(Float64, Float64)
    u = t - r
    v = t + r
    U = atan(u)
    V = atan(v)
    return (U,V)
end

function CompactifiedMinkowskiderivOP{T<:Int}(Nx::T, Ny::T, M::T, loc::Array{T,1})::Array{Float64, 4}
    @assert Nx == Ny
	operator = zeros(Nx+1, Ny+1, Nx+1, Ny+1)
    uglobal = chebgrid(Nx, M, loc[1])
    vglobal = chebgrid(Ny, M, loc[2])
    for index in CartesianRange(size(operator))    
        k = index.I[1]
        i = index.I[2]
        l = index.I[3]
        j = index.I[4]   
        ginvfac = ginv((pi/2)*uglobal[i], (pi/2)*vglobal[k])
        detgfac = detg((pi/2)*uglobal[i], (pi/2)*vglobal[k])
        operator[index] = ginvfac*detgfac*chebw(i,Ny)*chebw(k,Nx)*chebd(k,l,Nx)*chebd(i,j,Ny)	
	end
	return operator
end

function CompactifiedMinkowskicomputeonPatchBnd(fn::Function, Nx::Int, M::Int, loc::Int)::Array{Float64,1}
    fNpatch = Float64[fn(tan((pi/2)*U)) for U in chebgrid(Nx, M, loc)]
    return fNpatch  
end

function CompactifiedMinkowskigetPatchIC{T<:Integer}(fn::Function, s::T, Nx::T, M::T, loc::Int)::Boundary
    if s == 0
        return Boundary(0, CompactifiedMinkowskicomputeonPatchBnd(fn, Nx, M, loc))
    elseif s == 1
        return Boundary(1, CompactifiedMinkowskicomputeonPatchBnd(fn, Nx, M, loc))
    else
        error("Unknown direction passed.")
    end
end

function CompactifiedMinkowskidistribute{T<:Integer}(fbndr::Function, fbndc::Function, frhs::Function, Nx::T, Ny::T)::Dict{Array{Int,1}, Patch}
    M   = 2
    dop = derivOP(Nx, Ny)
    bop = boundaryOP(Nx, Ny)
    dbase  = Dict{Array{Int, 1}, Patch}()
    for i in 2:2M, k in i-min(i-1,M):min(i-1,M) 
        loc  = [k, i-k]
        if loc != [2,2]
            dop  = CompactifiedMinkowskiderivOP(Nx, Ny, M, loc)
            rhs  = RHS(frhs, Nx, Ny, M, loc)
            bndx = (loc[2]==1) ? (getPatchIC(fbndr, 0, Nx, M, loc[1])) : (getPatchBnd(dbase[loc-[0,1]], 0))
            bndy = (loc[1]==1) ? (getPatchIC(fbndc, 1, Ny, M, loc[2])) : (getPatchBnd(dbase[loc-[1,0]], 1))
            dbase[loc] = calcPatch(bndx, bndy, rhs, dop, bop, loc)
        else
            dbase[loc] = Patch(loc, zeros(Nx+1, Ny+1))
        end
    end
    return dbase
end
