#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 06-2018
#--------------------------------------------------------------------

function derivOP_corrected{T<:Int}(Nx::T, Ny::T)::Array{Float64, 4}
	" NOTE: This is the corrected version which evaluates the operator
			correctly when nx != ny. However, this doesn't seem to 
			give the correct solution."
	operator = zeros(Nx+1, Ny+1, Nx+1, Ny+1)
    for index in CartesianRange(size(operator)) 
        (m,n,p,q) = index.I
        operator[m,n,p,q] = 2*chebw(p,Nx)*chebw(n,Ny)*chebd(p,m,Nx)*chebd(n,q,Ny)
    end
    return operator
end

function derivOP{T<:Int}(Nx::T, Ny::T)::Array{Float64, 4}
    @assert Nx == Ny
	operator = zeros(Nx+1, Ny+1, Nx+1, Ny+1)
    for index in CartesianRange(size(operator))    
        k = index.I[1]
        i = index.I[2]
        l = index.I[3]
        j = index.I[4]   
        operator[index] = 2*chebw(i,Ny)*chebw(k,Nx)*chebd(k,l,Nx)*chebd(i,j,Ny)	
	end
	return operator
end

function boundaryOP{T<:Int}(Nx::T, Ny::T)::Array{Float64, 4}
    bnd = zeros(Nx+1, Ny+1, Nx+1, Ny+1)
    for index in CartesianRange(size(bnd))
        k = index.I[1]
        i = index.I[2]
        l = index.I[3]
        j = index.I[4]    
        if  i==1 || k==1
            bnd[index] = delta(i,j)*delta(k,l)
        end
    end
    return bnd
end

function RHS{T<:Int}(fn::Function, Nx::T, Ny::T, M::T, loc::Array{Int,1})::Array{Float64,2}
    rhs = projectonPatchbyRestriction(fn, Nx, Ny, M, loc)
    for index in CartesianRange(size(rhs))
        (i,j) = index.I
        rhs[i,j] = chebw(i,Nx)*chebw(j,Ny)*rhs[i,j]
    end
    return rhs
end
