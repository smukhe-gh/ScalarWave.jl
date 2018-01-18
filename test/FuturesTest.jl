#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------
# Functions to test a generic future implementation.
# Doen not test the functions in src 

@test typeof(@spawn 2) == Future

function fcomputeRHS(N::Int, M::Int, loc::Array{Int, 1}, fnbrow::Function, fnbcol::Function, dbase::Dict{Array{Int,1}, Future})::Future
     B2N  = initializeRHS(N, M, loc, fnbrow, fnbcol)
     if sum(loc) == 2
        return @spawn B2N
     elseif loc[1] == 1 || loc[2] == 1
        if loc[1] > loc[2]
            fB2N = fsetBC(B2N, fextractBC(dbase[loc-[1,0]], :R), :R)
        else
            fB2N = fsetBC(B2N, fextractBC(dbase[loc-[0,1]], :C), :C)
        end
     else
        fB2N = fsetBC(B2N, fextractBC(dbase[loc-[1,0]], :R), :R)
        fB2N = fsetBC(fetch(fB2N), fextractBC(dbase[loc-[0,1]], :C), :C)
     end
     return fB2N
end 

function pascal!(A::Array{Float64,2})::Array{Float64,2}
    N = size(A)[1]
    for i in 2:N, j in 2:N
        A[i,j] = A[i,j-1] + A[i-1,j]
    end
    return A
end

function fpascal(A::Future)::Future
    return @spawn pascal!(fetch(A))
end

function check_Futures(N::Int, M::Int)::Dict
    fdbase = Dict{Array{Int, 1}, Future}()
    for i in 2:2M, k in i-min(i-1,M):min(i-1,M)
        loc = [k, i-k]
        fRHS = fcomputeRHS(N-1, M, loc, (x,y)->1, (x,y)->1, fdbase)
        fdbase[loc] = fpascal(fRHS)
    end
    return fdbase
end

A  = ones(3,3)
fA = @spawn A
@test fetch(fpascal(fA)) == A

@test fetch(check_Futures(3,1)[[1,1]]) == pascal!(A)
