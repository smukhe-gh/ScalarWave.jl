#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

function guu(i, ii, j, jj, Px, Py, twist)
    return -cos((pi/twist)*cospi(chebx(i,Px)/2)*cospi(chebx(ii,Py)/2))*
            sin((pi/twist)*cospi(chebx(i,Px)/2)*cospi(chebx(ii,Py)/2))*
            delta(i,j)*delta(ii, jj)
end

function guv(i, ii, j, jj, Px, Py, twist)
    return +cos((pi/twist)*cospi(chebx(i,Px)/2)*cospi(chebx(ii,Py)/2))^2*delta(i,j)*delta(ii, jj)
end

function gvu(i, ii, j, jj, Px, Py, twist)
    return -sin((pi/twist)*cospi(chebx(i,Px)/2)*cospi(chebx(ii,Py)/2))^2*delta(i,j)*delta(ii, jj)
end

function gvv(i, ii, j, jj, Px, Py, twist)
    return +cos((pi/twist)*cospi(chebx(i,Px)/2)*cospi(chebx(ii,Py)/2))*
            sin((pi/twist)*cospi(chebx(i,Px)/2)*cospi(chebx(ii,Py)/2))*
            delta(i,j)*delta(ii, jj)
end

function W(i, ii, j, jj, Px, Py)
    return chebw(i, Px)*chebw(ii, Py)*delta(i,j)*delta(ii, jj) 
end

function dU(i, ii, j, jj, Px)
    return chebd(i, j, Px)*delta(ii, jj)  
end

function dV(i, ii, j, jj, Py)
    return delta(i, j)*chebd(ii, jj, Py)
end

function testderivOP(Px::Int, Py::Int)
    operator = zeros(Px+1, Py+1, Px+1, Py+1)
    arrW     = zeros(Px+1, Py+1, Px+1, Py+1)
    arrdU    = zeros(Px+1, Py+1, Px+1, Py+1)
    arrdV    = zeros(Px+1, Py+1, Px+1, Py+1)
    arrguu   = zeros(Px+1, Py+1, Px+1, Py+1)
    arrguv   = zeros(Px+1, Py+1, Px+1, Py+1)
    arrgvu   = zeros(Px+1, Py+1, Px+1, Py+1)
    arrgvv   = zeros(Px+1, Py+1, Px+1, Py+1)

    for index in CartesianRange(size(operator))
        l, ll, m, mm = index.I
        arrW[index]   = W(l, ll, m, mm, Px, Py) 
        arrdU[index]  = dU(l, ll, m, mm, Px) 
        arrdV[index]  = dV(l, ll, m, mm, Py) 
        arrguu[index] = 1.0 #guu(l, ll, m, mm, Px, Py, 1) 
        arrguv[index] = 1.0 #guv(l, ll, m, mm, Px, Py, 1) 
        arrgvu[index] = 1.0 #gvu(l, ll, m, mm, Px, Py, 1) 
        arrgvv[index] = 1.0 #gvv(l, ll, m, mm, Px, Py, 1) 
    end

    for index in CartesianRange(size(operator))
        l, ll, m, mm = index.I
        L1 = sum(arrW[i, ll, i, ll]*arrguu[i, ll, i, ll]*arrdU[i, ll, l, ll]*arrdU[i, ll, m, mm] 
                 for i in 1:Px+1) 
        L2 = sum(arrW[l, ii, l, ii]*arrgvv[l, ii, l, ii]*arrdV[l, ii, l, ll]*arrdV[l, ii, m, mm] 
                 for ii in 1:Py+1) 
        L3 = arrW[m, ll, m, ll]*arrguv[m, ll, m, ll]*arrdU[m, ll, l, ll]*arrdV[m, ll, m, mm] +
             arrW[l, mm, l, mm]*arrgvu[l, mm, l, mm]*arrdV[l, mm, l, ll]*arrdU[l, mm, m, mm]
        operator[index] = L3 # + L1 + L2
    end
    return operator
end

function testboundaryOP{T<:Int}(Px::T, Py::T)::Array{Float64, 4}
    bnd = zeros(Px+1, Py+1, Px+1, Py+1)
    for index in CartesianRange(size(bnd))
        i, ii, j, jj = index.I
        if  i==1 || ii==1
            bnd[index] = delta(ii,jj)*delta(i,j)
        end
    end
    return bnd
end

function testRHS(ufunc::Function, vfunc::Function, Px::Int, Py::Int)
    B = zeros(Px+1, Py+1)
    u = chebgrid(Px)
    v = chebgrid(Py)
    B[:,1] = ufunc.(u)
    B[1,:] = vfunc.(v)
    return B
end

Nx, Ny = 40, 40
@assert Nx == Ny
fH2L   = vec([x^3 + y^5 for x in chebgrid(Nx), y in chebgrid(Ny)])
opH2L  = shapeH2L(testderivOP(Nx, Ny))

# create a boundary operator that sets all the boundaries
#bpH2L  = shapeH2L(testboundaryOP(Nx, Ny))
Bmat   = zeros(Nx+1, Ny+1)
Bmat[1, :] = Bmat[end, :] = Bmat[:, end] = Bmat[:, 1] = 1
bpH2L  = diagm(vec(Bmat)) 
L      = opH2L + bpH2L

b      = shapeH2L(testRHS(u->0, v->10*exp(-v^2/0.2), Nx, Ny))

#= Operator construction
Du     = Float64[chebd(i, j, Nx) for i in 1:Nx+1, j in 1:Nx+1]
Dv     = Float64[chebd(i, j, Ny) for i in 1:Ny+1, j in 1:Ny+1]
Wu     = Float64[chebw(i, Nx) for i in 1:Nx+1]
Wv     = Float64[chebw(i, Ny) for i in 1:Ny+1]
WUV    = diagm(vec(kron(Wu, Wv)))
DV     = kron(eye(Nx+1), Dv)
DU     = kron(Du, eye(Ny+1))

twist  = 20
gUU    = diagm(vec(Float64[-cos((pi/twist)*cospi(chebx(i,Nx)/2)*cospi(chebx(j,Ny)/2))*
                            sin((pi/twist)*cospi(chebx(i,Nx)/2)*cospi(chebx(j,Ny)/2)) for i in 1:Nx+1, j in 1:Ny+1]))
gUV    = diagm(vec(Float64[+cos((pi/twist)*cospi(chebx(i,Nx)/2)*cospi(chebx(j,Ny)/2))^2 for i in 1:Nx+1, j in 1:Ny+1]))
gVU    = diagm(vec(Float64[-sin((pi/twist)*cospi(chebx(i,Nx)/2)*cospi(chebx(j,Ny)/2))^2 for i in 1:Nx+1, j in 1:Ny+1]))
gVV    = diagm(vec(Float64[+cos((pi/twist)*cospi(chebx(i,Nx)/2)*cospi(chebx(j,Ny)/2))*
                            sin((pi/twist)*cospi(chebx(i,Nx)/2)*cospi(chebx(j,Ny)/2)) for i in 1:Nx+1, j in 1:Ny+1])) 

L      = WUV*(gUU*DU*DU + gUV*DU*DV + gVU*DV*DU + gVV*DV*DV) 
=#


B      = bpH2L
u      = (L + B)\ b  
uL2H   = shapeL2H(u, Nx, Ny)

using PyPlot
plot(uL2H[end, :])
plot(uL2H[1, :])
#show()

@show cond(L + B)
@show sort(abs.(eigvals(L + B)))[1:3]

#drawmultipatch(Dict([1,1]=> Patch([1,1], uL2H)), "test-action-operator")
