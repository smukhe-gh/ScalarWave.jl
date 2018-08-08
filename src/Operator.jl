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


