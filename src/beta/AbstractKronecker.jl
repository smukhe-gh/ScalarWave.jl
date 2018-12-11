#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 11-2018
# Compute operator using Abstract Types
#--------------------------------------------------------------------

abstract type AbstractOperator{Tag} end

struct KroneckerProduct{Tag} <: AbstractOpeator{Tag}
    A::AbstractOperator
    B::AbstractOperator
end

