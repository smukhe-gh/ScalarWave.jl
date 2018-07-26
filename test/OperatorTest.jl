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
        arrguu[index] = guu(l, ll, m, mm, Px, Py, 1) 
        arrguv[index] = guv(l, ll, m, mm, Px, Py, 1) 
        arrgvu[index] = gvu(l, ll, m, mm, Px, Py, 1) 
        arrgvv[index] = gvv(l, ll, m, mm, Px, Py, 1) 
    end

    for index in CartesianRange(size(operator))
        l, ll, m, mm = index.I
        L1 = sum(arrW[i, ll, i, ll]*arrguu[i, ll, i, ll]*arrdU[i, ll, l, ll]*arrdU[i, ll, m, mm] 
                 for i in 1:Px+1) 
        L2 = sum(arrW[l, ii, l, ii]*arrgvv[l, ii, l, ii]*arrdV[l, ii, l, ll]*arrdV[l, ii, m, mm] 
                 for ii in 1:Py+1) 
        L3 = arrW[m, ll, m, ll]*arrguv[m, ll, m, ll]*arrdU[m, ll, l, ll]*arrdV[m, ll, m, mm] +
             arrW[l, mm, l, mm]*arrgvu[l, mm, l, mm]*arrdV[l, mm, l, ll]*arrdU[l, mm, m, mm]
        operator[index] = L1 + L2 + L3
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

fH2L  = vec([x^3 + y^5 for x in chebgrid(Nx), y in chebgrid(Ny)])
opH2L = shapeH2L(testderivOP(Nx, Ny), Nx, Ny)
bpH2L = shapeH2L(testboundaryOP(Nx, Ny), Nx, Ny)

# Now test the functions of both of these operators
@show fH2L'*opH2L*fH2L
@show bpH2L / fH2L
