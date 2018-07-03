#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 06-2018
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

#--------------------------------------------------------------------
# Now construct the derivative operators
#--------------------------------------------------------------------

# Consider the operator that works and the associated boundary operator
function derivOP{T<:Int}(Nx::T, Ny::T)::Array{Float64, 4}
    @assert Nx == Ny
	operator = zeros(Nx+1, Ny+1, Nx+1, Ny+1)
    for index in CartesianRange(size(operator))    
        (k,i,l,j) = index.I
        operator[index] = 2*chebw(i,Ny)*chebw(k,Nx)*chebd(k,l,Nx)*chebd(i,j,Ny)	
	end
	return operator
end

function boundaryOP{T<:Int}(Nx::T, Ny::T)::Array{Float64, 4}
    bnd = zeros(Nx+1, Ny+1, Nx+1, Ny+1)
    for index in CartesianRange(size(bnd))
        (k,i,l,j) = index.I
        if  i==1 || k==1
            bnd[index] = delta(i,j)*delta(k,l)
        end
    end
    return bnd
end

#=
Explicit operator computation should involve

    operator[l, ll, m, mm] = sum(quadW(i, ii, j, jj, Px, Py)*
                                (g00(i, ii, j, jj, Px, Py, twist)*derivU(i, ii, l, ll, Px)*derivU(i, ii, m, mm, Py) +  
                                 g01(i, ii, j, jj, Px, Py, twist)*derivU(i, ii, l, ll, Px)*derivV(i, ii, m, mm, Py) + 
                                 g10(i, ii, j, jj, Px, Py, twist)*derivV(i, ii, l, ll, Px)*derivU(i, ii, m, mm, Py) + 
                                 g11(i, ii, j, jj, Px, Py, twist)*derivV(i, ii, l, ll, Px)*derivV(i, ii, m, mm, Py)) 
                                 for i in 1:Px+1, ii in 1:Py+1, j in 1:Px+1, jj in 1:Py+1)

However, we should be able to simplify this due to the Kronecker deltas in each function
=#

# Now introduce the coordinate transformation, and compare with Mathematica
function derivOP_efficient(Px, Py, twist)
    operator = zeros(Px+1, Py+1, Px+1, Py+1)
    for index in CartesianRange(size(operator))
        l, ll, m, mm = index.I
        operator[l, ll, m, mm] = sum(quadW(i, ii, j, jj, Px, Py)*
                                (g00(i, ii, j, jj, Px, Py, twist)*derivU(i, ii, l, ll, Px)*derivU(i, ii, m, mm, Py) +  
                                 g01(i, ii, j, jj, Px, Py, twist)*derivU(i, ii, l, ll, Px)*derivV(i, ii, m, mm, Py) + 
                                 g10(i, ii, j, jj, Px, Py, twist)*derivV(i, ii, l, ll, Px)*derivU(i, ii, m, mm, Py) + 
                                 g11(i, ii, j, jj, Px, Py, twist)*derivV(i, ii, l, ll, Px)*derivV(i, ii, m, mm, Py)) 
                                 for i in 1:Px+1, ii in 1:Py+1, j in 1:Px+1, jj in 1:Py+1 if i==j && ii==jj)
    end
    return operator
end
# Compute the derivative operator with all summations included. This would 
# be incredibly slow for larger Px, Py
function derivop_copy(px, py)
    @assert px + py < 10 
    operator = zeros(px+1, py+1, px+1, py+1)
    for index in cartesianrange(size(operator))
        l, ll, m, mm = index.i
        operator[l, ll, m, mm] = sum( w(i, ii, j, jj, px, py)*
                                    (du(i, ii, l, ll, px)*du(i, ii, m, mm, py) +  
                                     du(i, ii, l, ll, px)*dv(i, ii, m, mm, py) + 
                                     dv(i, ii, l, ll, px)*du(i, ii, m, mm, py) + 
                                     dv(i, ii, l, ll, px)*dv(i, ii, m, mm, py)) 
                                     for i in 1:px+1, ii in 1:py+1, j in 1:px+1, jj in 1:py+1)
    end
    return operator
end

# Now introduce the coordinate transformation, and compare with Mathematica
function derivOP_efficient(Px, Py, twist)
    operator = zeros(Px+1, Py+1, Px+1, Py+1)
    for index in CartesianRange(size(operator))
        l, ll, m, mm = index.I
        operator[l, ll, m, mm] = sum(W(i, ii, i, ii, Px, Py)*
                                    (guu(i, ii, i, ii, Px, Py, twist)*dU(i, ii, l, ll, Px)*dU(i, ii, m, mm, Py) +  
                                     guv(i, ii, i, ii, Px, Py, twist)*dU(i, ii, l, ll, Px)*dV(i, ii, m, mm, Py) + 
                                     gvu(i, ii, i, ii, Px, Py, twist)*dV(i, ii, l, ll, Px)*dU(i, ii, m, mm, Py) + 
                                     gvv(i, ii, i, ii, Px, Py, twist)*dV(i, ii, l, ll, Px)*dV(i, ii, m, mm, Py)) 
                                     for i in 1:Px+1, ii in 1:Py+1)
    end
    return operator
end

# FIXME: Doesn't work after collapsing more indices
function derivOP_just_a_bit_more_efficient(Px, Py, twist)
    operator = zeros(Px+1, Py+1, Px+1, Py+1)
    for index in CartesianRange(size(operator))
        l, ll, m, mm = index.I
        operator[l, ll, m, mm] = sum(quadW(i, ll, i, ll, Px, Py)*g00(i, ll, i, ll, Px, Py, twist)*derivU(i, ll, l, ll, Px)*derivU(i, mm, m, mm, Py) 
                                     for i in 1:Px+1) + 
                                 sum(quadW(l, ii, l, ii, Px, Py)*g11(l, ii, l, ii, Px, Py, twist)*derivV(l, ii, l, ll, Px)*derivV(m, ii, m, mm, Py)
                                     for ii in 1:Py+1) +  
                                 quadW(m, ll, m, ll, Px, Py)*g01(m, ll, m, ll, Px, Py, twist)*derivU(m, ll, l, ll, Px)*derivV(m, ll, m, mm, Py) + 
                                 quadW(l, mm, l, mm, Px, Py)*g10(l, mm, l, mm, Px, Py, twist)*derivV(l, mm, l, ll, Px)*derivU(l, mm, m, mm, Py) 
    end
    return operator
    
end

#---------------------------------------------------------------------
# Do a coordinate transformation of the wave operator
#---------------------------------------------------------------------

function dalembertian_null_op(px, py)
    @assert px + py < 10 
    operator = zeros(px+1, py+1, px+1, py+1)
    for index in cartesianrange(size(operator))
        l, ll, m, mm = index.i
        operator[l, ll, m, mm] = sum((guv(i, ii, j, jj) + gvu(i,ii,j,jj))*
                                     dU(i, ii, j, jj, px)*dV(j, jj, m, mm, py) 
                                     for i in 1:px+1, ii in 1:py+1, j in 1:px+1, jj in 1:py+1)
    end
    return operator
end
