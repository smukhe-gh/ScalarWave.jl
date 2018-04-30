#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 04-2018
# Scalar Wave equation in fixed Schwarzschild background in 
# double null coordinates
#--------------------------------------------------------------------

function findrschwarzschild(u::Float64, v::Float64, MBH::Float64)::Float64
    f(r) = - exp(r/2MBH)*(r/2MBH -1) - u*v 
    return fzero(f, 2MBH)
end

function physcoords(umin::Float64, umax::Float64, Nx::Int, M::Int, loc::Int):Array{Float64,1}
    return Float64[(((umax - umin)/2)*x) + ((umax + umin)/2) for x in chebgrid(Nx, M, loc)]
end

function projectonPatchBndbyRestriction(fn::Function, umin::Float64, umax::Float64, Nx::Int, M::Int, loc::Int)::Array{Float64,1}
    f2Npatch = Float64[fn(x) for x in physcoords(umin, umax, 2Nx, M, loc)]
    invndmx  = inv(vandermonde(2Nx))
    pvndmx   = pseudovandermonde(2Nx, chebgrid(Nx))
    fNpatch  = pvndmx*restrictmodes!(invndmx*f2Npatch, Nx)
    return fNpatch  
end

function projectonPatchbyRestriction(fn::Function, umin::Float64, umax::Float64, vmin::Float64, vmax::Float64,  
                                     Nx::Int, Ny::Int, M::Int, loc::Array{Int,1})::Array{Float64,2}
    f2Npatch = Float64[fn(x,y) for x in physcoords(umin, umax, 2Nx, M, loc[1]), y in physcoords(vmin, vmax, 2Ny, M, loc[2])]
    invndmx  = inv(vandermonde(2Nx))
    invndmy  = inv(vandermonde(2Ny))
    pvndmx   = pseudovandermonde(2Nx, chebgrid(Nx))
    pvndmy   = pseudovandermonde(2Ny, chebgrid(Ny))
    fNpatch  = pvndmx*restrictmodes!(invndmx*f2Npatch*invndmy', Nx, Ny)*pvndmy'
    return fNpatch  
end

function getPatchIC{T<:Integer}(fn::Function, umin::Float64, umax::Float64, s::T, Nx::T, M::T, loc::Int)::Boundary
    if s == 0
        return Boundary(0, projectonPatchBndbyRestriction(fn, umin, umax, Nx, M, loc))
    elseif s == 1
        return Boundary(1, projectonPatchBndbyRestriction(fn, umin, umax, Nx, M, loc))
    else
        error("Unknown direction passed.")
    end
end

function detg(u::Float64, v::Float64, MBH::Float64)::Float64
    rsolve = findrschwarzschild(u, v, MBH)
    return ((32*MBH^3*exp(-rsolve/2MBH))/rsolve)^2
end

function derivOP{T<:Int}(umin::Float64, umax::Float64, vmin::Float64, vmax::Float64, MBH::Float64, 
                         Nx::T, Ny::T, M::T, loc::Array{Int,1})::Array{Float64, 4}
	operator = zeros(Nx+1, Ny+1, Nx+1, Ny+1)
    uspan    = physcoords(umin, umax, Nx, M, loc[1])
    vspan    = physcoords(vmin, vmax, Ny, M, loc[2])
    for index in CartesianRange(size(operator))    
        k = index.I[1]
        i = index.I[2]
        l = index.I[3]
        j = index.I[4]   
        operator[index] = 2*chebw(i,Ny)*chebw(k,Nx)*chebd(k,l,Nx)*chebd(i,j,Ny)*detg(uspan[i], uspan[k], MBH)^(3/2)
	end
	return operator
end

function RHS{T<:Int}(fn::Function, umin::Float64, umax::Float64, vmin::Float64, 
                     vmax::Float64, MBH::Float64,  Nx::T, Ny::T, M::T, loc::Array{Int,1})::Array{Float64,2}
    rhs = projectonPatchbyRestriction(fn, umin, umax, vmin, vmax, Nx, Ny, M, loc)
    uspan = physcoords(umin, umax, Nx, M, loc[1])
    vspan = physcoords(vmin, vmax, Ny, M, loc[2])
    for index in CartesianRange(size(rhs))
        i = index.I[1]
        j = index.I[2]
        rhs[i,j] = (chebw(i,Nx)/M)*(chebw(i,Ny)/M)*rhs[i,j]*detg(uspan[i], uspan[j], MBH)^(1/2)
    end
    return rhs
end

function distribute{T<:Integer}(fbndr::Function, fbndc::Function, frhs::Function, Nx::T, Ny::T, M::T, 
                                        umin::Float64, umax::Float64, vmin::Float64, vmax::Float64, MBH::Float64)::Dict{Array{Int,1}, Patch}
    bop = boundaryOP(Nx, Ny)
    u   = linspace(umax, umin, M+1)
    v   = linspace(vmax, vmin, M+1)
    dbase  = Dict{Array{Int, 1}, Patch}()
    for i in 2:2M, k in i-min(i-1,M):min(i-1,M) 
        loc  = [k, i-k]
        (lumin, lumax) = (u[k], u[k+1])
        (lvmin, lvmax) = (u[i-k], u[i-k+1])
        jacobianfac = ((lumax - lumin)*(lvmax - lvmin)/4)^2 
        dop  = derivOP(umin, lumax, lvmin, lvmax, MBH, Nx, Ny, M, loc)*jacobianfac
        rhs  = RHS(frhs, lumin, lumax, lvmin, lvmax, MBH, Nx, Ny, M, loc)
        @show loc, cond(shapeH2L(dop))
        bndx = (loc[2]==1) ? (getPatchIC(fbndr,  umin, umax, 0, Nx, M, loc[1])) : (getPatchBnd(dbase[loc-[0,1]], 0))
        bndy = (loc[1]==1) ? (getPatchIC(fbndc,  vmin, vmax, 1, Ny, M, loc[2])) : (getPatchBnd(dbase[loc-[1,0]], 1))
        dbase[loc] = calcPatch(bndx, bndy, rhs, dop, bop, loc)
    end
    return dbase
end
