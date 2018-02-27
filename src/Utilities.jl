#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

function delta{T<:Int}(i::T, j::T)::Float64
	return i==j ? 1 : 0
end

function coordtransL2G{T<:Int}(M::T, loc::T, xp::Float64)::Float64
    if loc > M || loc < 1
        error("Location incompatible with the number of patches")
    end
    s = Float64[d for d in M-1:-2:-M+1]
	xg = (xp + s[loc])/M
	return xg 
end

function coordtransG2L{T<:Int}(M::T, loc::T, xg::Float64)::Float64
    if loc > M || M < 1
        error("Location incompatible with the number of patches")
    end
    s = Float64[d for d in M-1:-2:-M+1]
	xp = xg*M - s[loc]
	return xp
end

function jacobian(M::Int)::Float64
    return 1/M
end

function shapeH2L(op4N::Array{Float64,4})::Array{Float64,2}
    Nx   = size(op4N)[1] - 1 
    Ny   = size(op4N)[2] - 1
    op2N = reshape(op4N, ((Nx+1)*(Ny+1),(Nx+1)*(Ny+1)))
    return op2N
end

function shapeH2L(b2N::Array{Float64,2})::Array{Float64,1}
    Nx = size(b2N)[1] - 1
    Ny = size(b2N)[2] - 1
    b1N = reshape(b2N, (Nx+1)*(Ny+1))
    return b1N
end

function shapeL2H(op2N::Array{Float64,2}, Nx::Int, Ny::Int)::Array{Float64,4}
    op4N = reshape(op2N, (Nx+1, Ny+1, Nx+1, Ny+1))
    return op4N
end

function shapeL2H(b1N::Array{Float64,1}, Nx::Int, Ny::Int)::Array{Float64,2}
    b2N = reshape(b1N, (Nx+1, Ny+1))
    return b2N
end

function LInfnorm{T<:Array{Float64,2}}(numericalGridData::T, exactGridData::T)::Float64
    errorGridData = numericalGridData - exactGridData
    return maximum(abs.(errorGridData))
end

function L1norm{T<:Array{Float64,2}}(numericalGridData::T, exactGridData::T, wx::Array{Float64,1}, wy::Array{Float64,1})::Float64
    errorGridData = numericalGridData - exactGridData
    return (wx'*(errorGridData)*wy)
end

function L2norm{T<:Array{Float64,2}}(numericalGridData::T, exactGridData::T, wx::Array{Float64,1}, wy::Array{Float64,1})::Float64
    errorGridData = numericalGridData - exactGridData
    return sqrt(wx'*(errorGridData.^2)*wy)
end

function dict2array(dbase::Dict{Array{Int,1}, Patch})::Array{Float64,2}
    Nx = size(dbase[[1,1]].value)[1] - 1
    Ny = size(dbase[[1,1]].value)[2] - 1
    M  = convert(Int, sqrt(length(dbase)))
    sPatch = zeros((Nx+1)*M, (Ny+1)*M)
    for m in 1:M, n in 1:M
        sPatch[1+(m-1)*(Nx+1):m*(Nx+1), 1+(n-1)*(Ny+1):n*(Ny+1)] = dbase[[m,n]].value
    end
    return sPatch
end

function array2dict(mPatch::Array{Float64,2}, Nx::Int, Ny::Int, M)::Dict{Array{Int,1}, Patch}
    dbase = Dict{Array{Int,1}, Patch}()
    for m in 1:M, n in 1:M
        dbase[[m,n]] = Patch([m,n], mPatch[1+(m-1)*(Nx+1):m*(Nx+1), 1+(n-1)*(Ny+1):n*(Ny+1)])
    end
    return dbase
end

function savegrid(dbase::Dict, path::String)
    datetime =  DateTime(now())
    jldopen("$path/$datetime.jld", "w") do file
        addrequire(file, ScalarWave)
        write(file, "patches", dbase)
    end
end

function loadgrid(path::String)::Dict
    dbase = load("/tmp/myfile.jld")
    return dbase
end

