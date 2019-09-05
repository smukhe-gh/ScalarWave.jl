#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

module ScalarWave

include("./Types/AbstractTypes.jl")
include("./Spectral/Basis/BasisTypes.jl")
include("./Spectral/Basis/ChebyshevGL.jl")
include("./Spectral/Spaces/1Dspace.jl")
include("./Spectral/Spaces/2Dspace.jl")
include("./Spectral/Spaces/AnySpace.jl")
include("./Visualization/PyPlot.jl")
include("./Misc/AxiSymmetry.jl")
include("./Misc/DoubleNullCoordinates.jl")
include("./Misc/NonLinSolverUtil.jl")
include("./Misc/InitialDataSolverUtil.jl")

end 
