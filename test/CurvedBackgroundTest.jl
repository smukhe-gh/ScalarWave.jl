#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

function testfindrschwarzschild(r::Float64, t::Float64, MBH::Float64)::Bool
    rt = r + 2MBH*log(abs(r/2MBH - 1))
    u  = t - rt
    v  = t + rt
    U  = r > 2MBH ?  -exp(-u/4MBH) : exp(-u/4MBH)
    V  = exp(v/4MBH)
    rsolve = findrschwarzschild(U, V, MBH)
    return rsolve ≈ r
end

function testdetg(r::Float64, t::Float64, MBH::Float64)::Bool
    rt = r + 2MBH*log(abs(r/2MBH - 1))
    u  = t - rt
    v  = t + rt
    U  = r > 2MBH ?  -exp(-u/4MBH) : exp(-u/4MBH)
    V  = exp(v/4MBH)
    return detg(U, V, MBH) ≈ ((32*MBH^3*exp(-r/2MBH))/r)^2
    
end

function testderivOP(umin::Float64, umax::Float64, vmin::Float64, vmax::Float64, MBH::Float64, Nx::Int, Ny::Int)::Array{Float64,2}
    wx = Float64[chebw(i,Nx) for i in 1:Nx+1]                   #size: (Nx+1)
    wy = Float64[chebw(i,Ny) for i in 1:Ny+1]                   #size: (Ny+1)
    dx = Float64[chebd(i,j,Nx) for i in 1:Nx+1, j in 1:Nx+1]    #size: (Nx+1, Nx+1)
    dy = Float64[chebd(i,j,Ny) for i in 1:Ny+1, j in 1:Ny+1]    #size: (Ny+1, Ny+1)
    
    ix = eye(Nx+1,Nx+1)     #size: (Nx+1, Nx+1)                                     
    iy = eye(Ny+1,Ny+1)     #size: (Ny+1, Ny+1)
    uspan = physcoords(umin, umax, Nx)
    vspan = physcoords(vmin, vmax, Ny)
    detgm = Float64[detg(u, v, MBH).^(3/2) for u in uspan, v in vspan] 
    detgmat = diagm(vec(detgm))
    w  = kron(wx,wy)        #size: (Nx+1, Ny+1)
    W  = diagm(vec(w))      #size: ((Nx+1)*(Ny+1), (Nx+1)*(Ny+1)) 
    DU = kron(iy, dx)       #size: ((Ny+1)*(Nx+1), (Ny+1)*(Nx+1)) 
    DV = kron(dy, ix)       #size: ((Ny+1)*(Nx+1), (Ny+1)*(Nx+1)) 
    WD = W*(DU*DV + DV*DU)  #size: ((Ny+1)*(Nx+1), (Ny+1)*(Nx+1)) 
    return detgmat*WD
end

@test testfindrschwarzschild(4.0, 3.0, 1.0) == 1
@test testfindrschwarzschild(1.0, 3.0, 1.0) == 1
@test testdetg(4.0, 3.45, 1.0) == 1

