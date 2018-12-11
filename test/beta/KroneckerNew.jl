#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 10-2018
# Test the structs for Kronecker products and the 
# the associated operations
#--------------------------------------------------------------------

import Base.getindex, Base. *, Base. +

abstract type AbstractOperator{S} end

struct IdentityOperator{S} <: AbstractOperator{S} end

struct KroneckerProduct{S1, S2, S} <: AbstractOperator{S}
    O1::AbstractOperator{S1}
    O2::AbstractOperator{S2}
end

# You should be able to 
#   - Multiply two Kronecker products
#   - Add to Kronecker products
#   - Compose higher order operations on the Kronecker Product
# All of these should be just representations, and no excess

