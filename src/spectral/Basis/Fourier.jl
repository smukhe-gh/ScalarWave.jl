#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 05-2019
# General Fourier Series: Endpoint Grid
# See Boyd F.2
#--------------------------------------------------------------------

function collocation(S::Fourier{Tag, N}, i::Int)::Real where {Tag, N}
    return (pi/N)*abs(i-1)
end

function basis(S::Fourier{Tag, N}, j::Int, x::Real)::Real where {Tag, N}
    return (1/2N)*sin(N*(x-collocation(S,j)))*cot(0.5*(x-collocation(S,j)))
end

function derivative(S::Fourier{Tag, N}, i::Int, j::Int)::Real where {Tag, N}
    if i != j
        return ((1/2)*((-1)^(i-j))*cot((1/2)*(collocation(S, i) - collocation(S, j))))
    else
        return 0
    end
end

function secondderivative(S::Fourier{Tag, N}, i::Int, j::Int)::Real where {Tag, N}
    if i != j
        return ((1/2)*((-1)^(i-j+1))/(sin((1/2)*(collocation(S, i) - collocation(S, j)))))
    else
        return -(1+2*N^2)/6
    end
end
