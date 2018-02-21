 #--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 02-2018
#--------------------------------------------------------------------

function rect(x,y,w,h)
    return Shape(x + [0,w,w,0], y + [0,0,h,h]) 
end
