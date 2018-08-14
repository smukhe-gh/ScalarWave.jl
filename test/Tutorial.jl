#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Tutorial 14.08.2018
#--------------------------------------------------------------------

# ==> Quick start
# To solve the wave equation on Minkowski (1+1) with null boundaries

patchdict = distribute(bnd1, bnd2, potential, Nx, Ny, M)

# where distribute (or fdistribute) if you want to use futures, 
# takes in the boundary conditions and any potential you might 
# want to specify, along with the grid arguments; i.e order of 
# polynomial approximation in each direction (Nx, Ny), and the
# number of patches in each direction. The function returns a 
# dictionary of patches, each with it's location (as the key) 
# and the solution in patch.value

# ==> The details
# If you want to solve a different equation you'd need to refine
# the function derivOP and/or the boundary operator. To see this, 
# let's look at the function body of distribute

function distribute{T<:Integer}(fbndr::Function, fbndc::Function,
                                frhs::Function, Nx::T, Ny::T, M::T)::Dict{Array{Int,1}, Patch}
    dop = derivOP(Nx, Ny)/M^4
    bop = boundaryOP(Nx, Ny)
    dbase  = Dict{Array{Int, 1}, Patch}()
    for i in 2:2M, k in i-min(i-1,M):min(i-1,M) 
        loc  = [k, i-k]
        rhs  = RHS(frhs, Nx, Ny, M, loc)
        bndx = (loc[2]==1) ? (getPatchIC(fbndr, 0, Nx, M, loc[1])) : (getPatchBnd(dbase[loc-[0,1]], 0))
        bndy = (loc[1]==1) ? (getPatchIC(fbndc, 1, Ny, M, loc[2])) : (getPatchBnd(dbase[loc-[1,0]], 1))
        dbase[loc] = calcPatch(bndx, bndy, rhs, dop, bop, loc)
    end
    return dbase
end

# While the boundary operator is simple 1's and 0's, 
# to compute the derivative operator you'd need lower 
# level functios that compute spectral differentiation matrices
# and integration weights. These are available in SpecCalculus.jl

function chebd{T<:Int}(i::T, j::T, N::T)::Float64
function chebw{T<:Int}(i::T, N::T)::Float64
function vandermonde(N::Int)::Array{Float64,2}

# You could change between the nodal representation of
# the solution to the modal representation
# or interpolate using the functions in Patch.jl

function extractPatchCoeffs(patch::Patch)::Array{Float64,2}
function interpolatePatch(patch::Patch, Nx::Int, Ny::Int)::Patch

# One also has the ability to refine or coarsen in h-p refinement,
# and the functions that let you do this are in Grid.jl

function restrictmodes!(coeffs::Array{Float64,2}, M::Int, N::Int)::Array{Float64,2}
function prolongatemodes(coeffs::Array{Float64,1}, M::Int)::Array{Float64,1}
function prolongatePatch(patch::Patch, M::Int)::Dict{Array{Int, 1}, Patch}
function restrictPatch(dbase::Dict{Array{Int, 1}, Patch})::Patch

# Finally, we have several functions in Utilities.jl that handle
# several auxillary functions such as array reshaping, computing norms etc

function coordtransL2G{T<:Int}(M::T, loc::T, xp::Float64)::Float64
function shapeH2L(op4N::Array{Float64,4})::Array{Float64,2}
function shapeH2L(b2N::Array{Float64,2})::Array{Float64,1}

# The script visualization.jl provides an interface to plot the solution
# using the Luxor.jl package.

drawmultipatch(dict, "visualization-test")
 
