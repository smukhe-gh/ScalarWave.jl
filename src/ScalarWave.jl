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
include("./Utilities/AxiSymmetry.jl")
include("./Utilities/DoubleNullCoordinates.jl")
include("./Utilities/Utilities.jl")
include("./Utilities/Physics.jl")
include("./Utilities/Diagnostics.jl")

end 
