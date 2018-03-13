#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 02-2018
#--------------------------------------------------------------------

function projectonPatchBndbyRestriction(fn::Function, Nx::Int)::Array{Float64,1}
    f2Npatch = Float64[fn(x) for x in chebgrid(2Nx)]
    invndmx  = inv(vandermonde(2Nx))
    pvndmx   = pseudovandermonde(2Nx, chebgrid(Nx))
    fNpatch  = pvndmx*restrictmodes!(invndmx*f2Npatch, Nx)
    return fNpatch  
end

function projectonPatchbyRestriction(fn::Function, Nx::Int, Ny::Int)::Array{Float64,2}
    f2Npatch = Float64[fn(x,y) for x in chebgrid(2Nx), y in chebgrid(2Ny)]
    invndmx  = inv(vandermonde(2Nx))
    invndmy  = inv(vandermonde(2Ny))
    pvndmx   = pseudovandermonde(2Nx, chebgrid(Nx))
    pvndmy   = pseudovandermonde(2Ny, chebgrid(Ny))
    fNpatch  = pvndmx*restrictmodes!(invndmx*f2Npatch*invndmy', Nx, Ny)*pvndmy'
    return fNpatch  
end

function projectonPatchBndbyRestriction(fn::Function, Nx::Int, M::Int, loc::Int)::Array{Float64,1}
    f2Npatch = Float64[fn(x) for x in chebgrid(2Nx, M, loc)]
    invndmx  = inv(vandermonde(2Nx))
    pvndmx   = pseudovandermonde(2Nx, chebgrid(Nx))
    fNpatch  = pvndmx*restrictmodes!(invndmx*f2Npatch, Nx)
    return fNpatch  
end

#=
#FIXME: We are doing something inconsistent here. 
function projectonPatchbyRestriction(fn::Function, Nx::Int, Ny::Int, M::Int, loc::Array{Int,1}, nd::Int, pd::Int)::Array{Float64,2}
    N2x, N2y = (2Nx+1, 2Ny+1)
    f2Npatch = Float64[fn(x,y) for x in chebgrid(N2x, M, loc[1]), y in chebgrid(N2y, M, loc[2])]
    invndmx  = inv(vandermonde(2Nx+1))
    invndmy  = inv(vandermonde(2Ny+1))
    fmodal   = invndmx*f2Npatch*invndmy'
    vndmx    = Float64[cheb(m, x) for x in chebgrid(Nx+nd), m in 0:Nx+pd]
    vndmy    = Float64[cheb(m, x) for x in chebgrid(Ny+nd), m in 0:Ny+pd]
    rfmodal  = fmodal[1:Nx+1+pd, 1:Ny+1+pd]
    fNpatch  = vndmx*rfmodal*vndmy'
    return fNpatch  
end
=#

function projectonPatchbyRestriction(fn::Function, Nx::Int, Ny::Int, M::Int, loc::Array{Int,1})::Array{Float64,2}
    f2Npatch = Float64[fn(x,y) for x in chebgrid(2Nx, M, loc[1]), y in chebgrid(2Ny, M, loc[2])]
    invndmx  = inv(vandermonde(2Nx))
    invndmy  = inv(vandermonde(2Ny))
    pvndmx   = pseudovandermonde(2Nx, chebgrid(Nx))
    pvndmy   = pseudovandermonde(2Ny, chebgrid(Ny))
    fNpatch  = pvndmx*restrictmodes!(invndmx*f2Npatch*invndmy', Nx, Ny)*pvndmy'
    return fNpatch  
end
