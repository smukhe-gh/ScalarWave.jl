#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

# can an array be assembled with futures?
function fextractBC(fpatch::Any, s::Int, N::Int):
	s == 0 ? @spawn extractBC(fetch(fpatch)), 0, N) : @spawn extractBC(fetch(fpatch)), 1, N)
end

function fsetBC(fpatch::Any, fboundary::Any, s::Int, N::Int)
	s == 0 ? @spawn setBC!(fetch(fpatch), fetch(fboundary), 0, N) : @spawn setBC!(fetch(fpatch), fetch(fboundary), 1, N)
end

