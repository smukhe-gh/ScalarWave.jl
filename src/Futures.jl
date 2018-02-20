#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

# TODO: These are not the functions we wrap. See notes. 
function fextractBC(patch::Future, s::Symbol)::Future
    @spawn extractBC(fetch(patch), s)
end

function fsetBC(patch::Array{Float64,2}, patchBC::Future, s::Symbol)::Future
    @spawn setBC(patch, fetch(patchBC), s)
end
 
