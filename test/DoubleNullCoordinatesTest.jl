#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
# Coordinate transformation Test
#--------------------------------------------------------------------

@testset "coordinate transformation" begin
    for i in 1:3, j in 1:3
        u = -abs(4*randn())
        v =  abs(4*randn())

        t = find_t_of_UV(u,v,1.0)
        r = find_r_of_UV(u,v,1.0)
        u_of_tr = find_U_of_tr(t,r,1.0) 
        v_of_tr = find_V_of_tr(t,r,1.0) 

        @test u ≈ u_of_tr
        @test v ≈ v_of_tr
    end
end

#--------------------------------------------------------------------
# Compare the coordinates with Mathematica
#--------------------------------------------------------------------
M, ω = (1.0, 1.0)
PV, PU = 100, 100
Umax, Umin = -4M, -8M
Vmin, Vmax =  4M,  8M
SUV = ProductSpace{GaussLobatto(V,PV, Vmax, Vmin), 
                   GaussLobatto(U,PU, Umax, Umin)}

#--------------------------------------------------------------------
# Define coordinates and their associated derivatives
#--------------------------------------------------------------------
𝕌 = Field(SUV, (U,V)->U)
𝕍 = Field(SUV, (U,V)->V)
θ = Field(SUV, (U,V)->π/2)
ϕ = Field(SUV, (U,V)->0)

𝕥  = Field(SUV, (U,V)->find_t_of_UV(U, V, M), 𝕌, 𝕍)
𝕣  = Field(SUV, (U,V)->find_r_of_UV(U, V, M), 𝕌, 𝕍)
NU = Field(SUV, (t,r)->find_U_of_tr(t, r, M), 𝕥, 𝕣)
NV = Field(SUV, (t,r)->find_V_of_tr(t, r, M), 𝕥, 𝕣)

@show maximum(abs(NU - 𝕌))
@show maximum(abs(NV - 𝕍))

drawpatch(abs(NU-𝕌), "coordinates-U")
drawpatch(abs(NV-𝕍), "coordinates-V")

using HDF5

if isfile("../output/hdf5/coordinate-transformation-collocation-points.h5")
    println("File already exits. Skipping")
else
    println("Creating dataset.")
    h5open("../output/hdf5/coordinate-transformation-collocation-points.h5", "w") do file
        write(file, "collocation-points-U",  𝕌.value)
        write(file, "collocation-points-V",  𝕍.value)
    end
end

# Read t and r from Mathematica
if isfile("../output/hdf5/coordinate-values-for-julia.h5")
    mathT = Field(SUV, h5read("../output/hdf5/coordinate-values-for-julia.h5", "t"))
    mathR = Field(SUV, h5read("../output/hdf5/coordinate-values-for-julia.h5", "r"))
else
    println("Could not find file. Create them using Mathematica")
    exit()
end

@show maximum(abs(mathT - 𝕥))
@show maximum(abs(mathR - 𝕣))

drawpatch(abs(mathT - 𝕥), "coordinate-error-t")
drawpatch(abs(mathR - 𝕣), "coordinate-error-r")

