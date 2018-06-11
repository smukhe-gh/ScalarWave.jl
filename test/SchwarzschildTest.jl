#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 05-2018
#--------------------------------------------------------------------

parameters = Dict("px"   => 20,
                  "py"   => 20,
                  "umin" => -20.0,
                  "umax" => -0.1,
                  "vmin" => 0.1,
                  "vmax" => 20.0,
                  "mass" => 1.0)

params  = dict2struct(parameters)
grid    = setgrid(params)
varlist = setvarlist(grid) 
dXdU    = 2/(grid.params.dmax[1] - grid.params.dmin[1])
dXdV    = 2/(grid.params.dmax[2] - grid.params.dmin[2])

# test coordinates
for index in CartesianRange(params.size)
    U = grid.U[index]
    V = grid.V[index]
    t = grid.t[index]
    r = grid.r[index]
    @test U ≈ find_UV_of_TR(t,r, grid.params.mass)[1] 
    @test V ≈ find_UV_of_TR(t,r, grid.params.mass)[2] 
end

# test incoming scalarwave
dict = SchwarzschildDistribute(u->0, v->exp(-v^2/0.1), (u,v)->0, parameters)

# plot the solution
drawmultipatch(dict, "schwarzschild-test-p-15")

# look at the eigen values and condition number of the operator
dop = SchwarzschildDerivOP(grid, varlist)
dop_old      = derivOP(20, 20)
bop          = boundaryOP(grid.params.mode[1], grid.params.mode[2])
operator     = shapeH2L(dop + bop)
operator_old = shapeH2L(dop_old + bop)
@show cond(operator)
@show cond(operator_old)

# Do a p-refinement and check for self-convergence
println("--------------------------------------------")
println(" Doing a self-convergence test")
println("--------------------------------------------")
for nmodes in 2:2:20
    parameters["px"] = parameters["py"] = nmodes
    params  = dict2struct(parameters)
    grid    = setgrid(params)
    varlist = setvarlist(grid) 
    dict1  = SchwarzschildDistribute(u->0, v->exp(-v^2/0.1), (u,v)->0, parameters)
    patch1 = interpolatePatch(dict1[[1,1]], 4*nmodes, 4*nmodes).value 

    parameters["px"] = parameters["py"] = nmodes + 1
    params  = dict2struct(parameters)
    grid    = setgrid(params)
    varlist = setvarlist(grid) 
    dict2  = SchwarzschildDistribute(u->0, v->exp(-v^2/0.1), (u,v)->0, parameters)
    patch2 = interpolatePatch(dict2[[1,1]], 4*nmodes, 4*nmodes).value 

    L2 = L2norm(patch1, patch2, chebweights(4*nmodes)*dXdU, chebweights(4*nmodes)*dXdV)^2
    @printf("p = %s L2 = %e \n", lpad(nmodes,3," "), L2) 
end

# Do a p-refinement and check for convergence
# with the wrong analytic solution
println("--------------------------------------------")
println(" Doing a fake-convergence test")
println("--------------------------------------------")
for nmodes in 2:2:10
    parameters["px"] = parameters["py"] = nmodes
    params  = dict2struct(parameters)
    grid    = setgrid(params)
    varlist = setvarlist(grid) 
    dict1  = SchwarzschildDistribute(u->0, v->exp(-v^2/0.1), (u,v)->0, parameters)
    patch1 = interpolatePatch(dict1[[1,1]], 4*nmodes, 4*nmodes).value 

    (Nx, Ny, M) = (nmodes, nmodes, 1)
    dict2  = distribute(u->0, v->exp(-v^2/0.1), (u,v)->0, Nx, Ny, M)
    patch2 = interpolatePatch(dict2[[1,1]], 4*nmodes, 4*nmodes).value 

    L2 = L2norm(patch1, patch2, chebweights(4*nmodes)*dXdU, chebweights(4*nmodes)*dXdV)^2
    @printf("p = %s L2 = %e \n", lpad(nmodes,3," "), L2) 
end
