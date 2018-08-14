#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

# FIXME : Check the metric functions
function fguu{T<:Int}(i::T, j::T, N::T, pibytwist::Float64)::Float64
    return -cos((pibytwist)*cospi(chebx(i,N)/2)*cospi(chebx(j,N)/2))*
            sin((pibytwist)*cospi(chebx(i,N)/2)*cospi(chebx(j,N)/2))
end

function fguv{T<:Int}(i::T, j::T, N::T, pibytwist::Float64)::Float64
    return +cos((pibytwist)*cospi(chebx(i,N)/2)*cospi(chebx(j,N)/2))*
            cos((pibytwist)*cospi(chebx(i,N)/2)*cospi(chebx(j,N)/2))
end

function fgvu{T<:Int}(i::T, j::T, N::T, pibytwist::Float64)::Float64
    return -cos((pibytwist)*cospi(chebx(i,N)/2)*cospi(chebx(j,N)/2))*
            cos((pibytwist)*cospi(chebx(i,N)/2)*cospi(chebx(j,N)/2))
end

function fgvv{T<:Int}(i::T, j::T, N::T, pibytwist::Float64)::Float64
    return +cos((pibytwist)*cospi(chebx(i,N)/2)*cospi(chebx(j,N)/2))*
            sin((pibytwist)*cospi(chebx(i,N)/2)*cospi(chebx(j,N)/2))
end

N   = 2
pibytwist = 0.0
Du  = Float64[chebd(i, j, N) for i in 1:N+1, j in 1:N+1]
Dv  = Float64[chebd(i, j, N) for i in 1:N+1, j in 1:N+1]
Wu  = Float64[chebw(i, N) for i in 1:N+1]
guu = Float64[fguu(i,j,N, pibytwist) for i in 1:N+1, j in 1:N+1]
guv = Float64[fguv(i,j,N, pibytwist) for i in 1:N+1, j in 1:N+1]
gvu = Float64[fgvu(i,j,N, pibytwist) for i in 1:N+1, j in 1:N+1]
gvv = Float64[fgvv(i,j,N, pibytwist) for i in 1:N+1, j in 1:N+1]
Wv  = Float64[chebw(i, N) for i in 1:N+1]

ginvuu = zeros(N+1, N+1)
ginvuv = zeros(N+1, N+1)
ginvvu = zeros(N+1, N+1)
ginvvv = zeros(N+1, N+1)
detg   = zeros(N+1, N+1)

for i in 1:N+1, j in 1:N+1 
    gmatrix = [guu[i,j] guv[i,j] 0  0;
              -gvu[i,j] gvv[i,j] 0  0;
               0        0        1  0;
               0        0        0  1]
    @show gmatrix
    @show ginv
    ginvuu[i,j] = ginv[1,1]
    ginvuv[i,j] = ginv[1,2]
    ginvvu[i,j] = ginv[2,1]
    ginvvv[i,j] = ginv[2,2]
    detg[i,j]   = det(gmatrix)
end

# Compute auxillary quantities
WUV = diagm(vec(kron(Wu, Wv)))
DV  = kron(eye(N+1), Dv)
DU  = kron(Du, eye(N+1))

ginvUU = diagm(vec(ginvuu))
ginvUV = diagm(vec(ginvuv))
ginvVU = diagm(vec(ginvvu))
ginvVV = diagm(vec(ginvvv))
detG   = diagm(vec(detg))

dginvUUU = diagm(DU*vec(ginvuu))
dginvUUV = diagm(DV*vec(ginvuu))
dginvUVU = diagm(DU*vec(ginvuv))
dginvUVV = diagm(DV*vec(ginvuv))
dginvVUU = diagm(DU*vec(ginvvu))
dginvVUV = diagm(DV*vec(ginvvu))
dginvVVU = diagm(DU*vec(ginvvv))
dginvVVV = diagm(DV*vec(ginvvv))
ddetgU   = diagm(DU*vec(detg))
ddetgV   = diagm(DV*vec(detg))

# Operator computation with derivatives of the determinant
L1 = ginvUU*DU*DU + ginvUV*DU*DV + ginvVU*DV*DU + ginvVV*DV*DV
   + dginvUUU*DU + dginvUVU*DV + dginvVUV*DU + dginvVVV*DV  
   + (1./2.(detG))*(ddetgU*ginvUU*DU + ddetgU*ginvUV*DV + ddetgV*ginvVU*DU + ddetgV*ginvVV*DV) 


# Compute the operator [long version] (L1 == L2)
L2 = ginvUU*DU*DU + ginvUV*DU*DV + ginvVU*DV*DU + ginvVV*DV*DV
+ dginvUUU*DU + (1/2)*dginvUUU*ginvUU*ginvUU*DU + (1/2)*dginvUUU*ginvUU*ginvUV*DV 
+ dginvUVU*DV + (1/2)*dginvUVU*ginvUU*ginvUV*DU + (1/2)*dginvUVU*ginvUV*ginvUV*DV
              + (1/2)*dginvVUU*ginvUU*ginvVU*DU + (1/2)*dginvVUU*ginvUV*ginvVU*DV
              + (1/2)*dginvVVU*ginvUU*ginvVV*DU + (1/2)*dginvVVU*ginvUV*ginvVV*DV
              + (1/2)*dginvUUV*ginvUU*ginvVU*DU + (1/2)*dginvUUV*ginvUU*ginvVV*DV
              + (1/2)*dginvUVV*ginvUV*ginvVU*DU + (1/2)*dginvUVV*ginvUV*ginvVV*DV
+ dginvVUV*DU + (1/2)*dginvVUV*ginvVU*ginvVU*DU + (1/2)*dginvVUV*ginvVU*ginvVV*DV
+ dginvVVV*DV + (1/2)*dginvVVV*ginvVU*ginvVV*DU + (1/2)*dginvVVV*ginvVV*ginvVV*DV

# Now solve the system
M     = 1
Nx    = Ny = N
dop   = shapeL2H(L1, N, N)/M^4
bop   = boundaryOP(Nx, Ny)
dbase = Dict{Array{Int, 1}, Patch}()
frhs  = (u,v) -> 0 
fbndr =  u    -> exp(-u^2/0.1) 
fbndc =  v    -> exp(-v^2/0.1) 

for i in 2:2M, k in i-min(i-1,M):min(i-1,M) 
    loc  = [k, i-k]
    rhs  = RHS(frhs, Nx, Ny, M, loc)
    bndx = (loc[2]==1) ? (getPatchIC(fbndr, 0, Nx, M, loc[1])) : (getPatchBnd(dbase[loc-[0,1]], 0))
    bndy = (loc[1]==1) ? (getPatchIC(fbndc, 1, Ny, M, loc[2])) : (getPatchBnd(dbase[loc-[1,0]], 1))
    dbase[loc] = calcPatch(bndx, bndy, rhs, dop, bop, loc)
end

drawmultipatch(dbase, "test-distorted-minkowski")
