#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

function delta(i::T, j::T)::Float64 where {T}
	return i==j ? 1 : 0
end

function LinearAlgebra. norm(f::Field{S})::Float64 where {S}
    return norm(reshape(f));
end

function writevtk(u::Field{ProductSpace{S1, S2, S3}}, filename::String) where {S1, S2, S3}
    length(S1) == length(S2) == length(S3) ? 0 : @warn "Paraview can act strange if the number of points are not the same in all directions"
    C1 = Field(S1, x->x)
    C2 = Field(S2, x->x)
    C3 = Field(S3, x->x)
    vtkfile = vtk_grid(filename, C1.value, C2.value, C3.value) 
    vtk_point_data(vtkfile, u.value, "scalar")
    outfiles = vtk_save(vtkfile)
end

function writevtk(u::Field{ProductSpace{S1, S2, S3}}, filename::String) where {S1, S2, S3 <: Chebyshev{Tag}} where {Tag}
    # FIXME: Hack; clean this up. 
    @assert length(S1) == length(S2) == length(S3)
    C1 = Float64.(Field(Taylor{U, order(S1), length(S1), 0}, x->x).value)
    vtkfile = vtk_grid(filename, C1, C1, C1) 
    vtk_point_data(vtkfile, u.value, "scalar")
    outfiles = vtk_save(vtkfile)
end
