#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

function derivOP{T<:Int}(Nx::T, Ny::T)::Array{Float64, 4}
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

function RHS{T<:Int}(Nx::T, Ny::T, fn::Function)::Array{Float64,2}
    return zeros(Nx+1, Ny+1) 
end

