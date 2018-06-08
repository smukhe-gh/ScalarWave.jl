#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 04-2018
# Compactified coordinates Metric functions
#--------------------------------------------------------------------

function power(x, a)   
    return x^a
end

function dict2struct(params::Dict)::Params
    return Params((params["px"], params["py"]),
                  (params["px"] + 1, params["py"] + 1),
                  (params["umin"], params["vmin"]),
                  (params["umax"], params["vmax"]),
                   params["mass"])
end

function find_TR_of_UV(U::Float64, V::Float64)::Tuple
    @vars x
    u = tan(find_zero(atan((tan(x)/(sqrt(1+tan(x)^2)))*log(1+tan(x)^2)) - U, (-pi/2, pi/2)))
    v = tan(find_zero(atan((tan(x)/(sqrt(1+tan(x)^2)))*log(1+tan(x)^2)) - V, (-pi/2, pi/2)))
    t = u+v
    r = v-u
    return (t, r)
end

function find_UV_of_TR(t::Float64, r::Float64)::Tuple
    u = (t-r)/2
    v = (t+r)/2
    U = atan((u/(sqrt(1+u^2)))*log(1+u^2))
    V = atan((v/(sqrt(1+v^2)))*log(1+v^2))
    return (U,V)
end

function dX_of_var(var::Array{Float64,2}, grid::Grid, X::Direction)
    dvar = zeros(grid.params.size)
    dXdU = 2/(grid.params.xmax[Int(X)] - grid.params.xmin[Int(X)])
    for index in CartesianRange(grid.params.size)
        (i, j) = index.I
        if Int(X) == 1
            dvar[index] = dXdU*sum(chebd(i, m, grid.params.p[Int(X)])*var[m, j] for m in 1:grid.params.size[Int(X)])
        else
            dvar[index] = dXdU*sum(chebd(j, n, grid.params.p[Int(X)])*var[i, n] for n in 1:grid.params.size[Int(X)])
        end
    end
    return dvar
end

function ddX_of_var(var::Array{Float64,2}, grid::Grid, X1::Direction, X2::Direction)
    ddvar = zeros(grid.params.size)
    dXdU  = 2/(grid.params.xmax[Int(X1)] - grid.params.xmin[Int(X1)])
    dXdV  = 2/(grid.params.xmax[Int(X2)] - grid.params.xmin[Int(X2)])
    for index in CartesianRange(grid.params.size)
        i, j = index[1], index[2]
        if (Int(X1) == 1 && Int(X2) == 1)
            ddvar[index] = dXdU*dXdU*sum(chebd(i, l, grid.params.p[Int(X1)])*chebd(l, m, grid.params.p[Int(X1)])*var[m,j] 
                                         for m in 1:grid.params.size[Int(X1)], l in 1:grid.params.size[Int(X1)])
        elseif (Int(X1) == 2 && Int(X2) == 2)
            ddvar[index] = dXdV*dXdV*sum(chebd(j, k, grid.params.p[Int(X2)])*chebd(k, n, grid.params.p[Int(X2)])*var[i,n] 
                                         for n in 1:grid.params.size[Int(X2)], k in 1:grid.params.size[Int(X2)])
        else
            ddvar[index] = dXdU*dXdV*sum(chebd(j, n, grid.params.p[Int(X1)])*chebd(i, k, grid.params.p[Int(X2)])*var[k,n] 
                                         for k in 1:grid.params.size[Int(X1)], n in 1:grid.params.size[Int(X2)])
        end
    end
    return ddvar
end
