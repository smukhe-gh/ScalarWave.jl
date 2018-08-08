#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

function fgetPatchIC{T<:Integer}(fbndr::Function, s::T, N::T, M::T, loc::T)::Future
    @spawnat workers()[rand(1:end)] getPatchIC(fbndr, s, N, M, loc)
end

function fRHS{T<:Integer}(frhs::Function, Nx::T, Ny::T, M::T, loc::Array{Int,1})::Future
    @spawnat workers()[rand(1:end)] RHS(frhs, Nx, Ny, M, loc)
end

function fgetPatchBnd(fpatch::Future, s::Int)::Future
    @spawnat workers()[rand(1:end)] getPatchBnd(fetch(fpatch), s)
end

function fcalcPatch{T<:Array{Float64,4}}(fbndx::Future, fbndy::Future, frhs::Future, dop::T, bop::T, loc::Array{Int,1})::Future
    @spawnat workers()[rand(1:end)] calcPatch(fetch(fbndx), fetch(fbndy), fetch(frhs), dop, bop, loc)
end
