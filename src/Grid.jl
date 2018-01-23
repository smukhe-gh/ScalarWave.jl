#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

function distribute{T<:Integer}(N::T, M::T)::Dict
	op = operator(N, M)
	dbase = Dict{Array{Int, 1}, Patch}()
    for i in 2:2M, k in i-min(i-1,M):min(i-1,M) 
        loc  = [k, i-k]	
        brow = loc[1]==1 ? getIC(loc, :R) : getPB(dbase[loc-[1,0], :R)
        bcol = loc[2]==1 ? getIC(loc, :C) : getPB(dbase[loc-[0,1], :C)
        boundary   = Boundary(brow, bcol)
        dbase[loc] = calcPatch(loc, boundary, op)
	end	
	return dbase
end

function fdist(N)
    patch = Array{Future}(N+1,N+1)
    for j in 1:N+1, i in 1:N+1
        fbnd0 = if i==1 fgetIC(0, i,j) else fgetPB(0, fpatch[i-1,j]) end
        fbnd1 = ...
        fpatch[i,j] = fcalcP(i,j, fbnd0,bnd1)
    end
    fpatch
end

fgetIC(d, i,j) = @spawn getIC(d, i,j)
fgetPB(d, fpatch) = @spawn getPB(d, fetch(fpatch))
fcalcP(i,j, fbnd0,bnd1) = @spawn calcP(i,j, fetch(bnd0),fetch(bnd1))

