#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 05-2018
# Type structure to handle metrics
#--------------------------------------------------------------------

struct Metric{T<:Union{Array{Float64,4}, Float64}}
    g00::T
    g11::T
    g22::T
    g33::T
end

function convertSoA2AoS(mesh::Mesh, SoAmetric::Metric)::Array{ScalarWave.Metric,4}
    AoSmetric = Array{Metric}(mesh.n) 
    for index in CartesianRange(mesh.n)
        AoSmetric[index] = Metric(SoAmetric.g00[index], SoAmetric.g11[index], 
                                  SoAmetric.g22[index], SoAmetric.g33[index]) 
    end
    return AoSmetric
end

function convertAoS2SoA(mesh::Mesh, AoSmetric::Array{ScalarWave.Metric,4})::Metric 
    SoAmetric = Metric(Array{Float64,4}(mesh.n), 
                       Array{Float64,4}(mesh.n), 
                       Array{Float64,4}(mesh.n), 
                       Array{Float64,4}(mesh.n))
    for index in CartesianRange(mesh.n)
        SoAmetric.g00[index] = AoSmetric[index].g00 
        SoAmetric.g11[index] = AoSmetric[index].g11 
        SoAmetric.g22[index] = AoSmetric[index].g22 
        SoAmetric.g22[index] = AoSmetric[index].g33 
    end
    return SoAmetric
end

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
