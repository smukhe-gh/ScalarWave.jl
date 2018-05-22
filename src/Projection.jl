#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 02-2018
#--------------------------------------------------------------------

function projectonPatchBndbyRestriction(fn::Function, Nx::Int, M::Int, loc::Int)::Array{Float64,1}
    f2Npatch = Float64[fn(x) for x in chebgrid(2Nx, M, loc)]
    invndmx  = inv(vandermonde(2Nx))
    pvndmx   = pseudovandermonde(2Nx, chebgrid(Nx))
    fNpatch  = pvndmx*restrictmodes!(invndmx*f2Npatch, Nx)
    return fNpatch  
end

function projectonPatchbyRestriction(fn::Function, Nx::Int, Ny::Int, M::Int, loc::Array{Int,1})::Array{Float64,2}
    f2Npatch = Float64[fn(x,y) for x in chebgrid(2Nx, M, loc[1]), y in chebgrid(2Ny, M, loc[2])]
    invndmx  = inv(vandermonde(2Nx))
    invndmy  = inv(vandermonde(2Ny))
    pvndmx   = pseudovandermonde(2Nx, chebgrid(Nx))
    pvndmy   = pseudovandermonde(2Ny, chebgrid(Ny))
    fNpatch  = pvndmx*restrictmodes!(invndmx*f2Npatch*invndmy', Nx, Ny)*pvndmy'
    return fNpatch  
end
