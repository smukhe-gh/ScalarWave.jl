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
    return (1/M)^2
end

function reshapeA(op4N::Array{Float64,4})::Array{Float64,2}
    N = size(op4N)[1] - 1
    op2N = reshape(op4N, ((N+1)^2, (N+1)^2))
    return op2N
end

function shapeA(op2N::Array{Float64,2})::Array{Float64,4}
    N = Int(sqrt(size(op2N)[1])) - 1
    op4N = reshape(op2N, (N+1, N+1, N+1, N+1))
    return op4N
end

function reshapeB(b2N::Array{Float64,2})::Array{Float64,1}
    N = size(b2N)[1] - 1
    b1N = reshape(b2N, (N+1)^2)
    return b1N
end

function shapeB(b1N::Array{Float64,1})::Array{Float64,2}
    N = Int(sqrt(size(b1N)[1])) - 1
    b2N = reshape(b1N, (N+1, N+1))
    return b2N
end

function LInfnorm{T<:Array{Float64,2}}(numericalGridData::T, exactGridData::T)::Float64
    errorGridData = numericalGridData - exactGridData
    return maximum(abs.(errorGridData))
end

function L1norm{T<:Array{Float64,2}}(numericalGridData::T, exactGridData::T, w::Array{Float64,1})::Float64
    errorGridData = numericalGridData - exactGridData
    return (w'*(errorGridData)*w)/(w'*(exactGridData)*w)
end

function L2norm{T<:Array{Float64,2}}(numericalGridData::T, exactGridData::T, w::Array{Float64,1})::Float64
    errorGridData = numericalGridData - exactGridData
    return sqrt((w'*(errorGridData.^2)*w)/(w'*(exactGridData.^2)*w))
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

