#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2018
#--------------------------------------------------------------------

"""
We want to be able to define two operations at 
the moment.
    -- Multiplication of two fields [Done]
    -- Derivative of a field [Done]
    -- Go between spaces
"""

struct M end
struct MU end
struct MV end
struct Mθ end

S  = GaussLobatto{M}(9)
x  = Field(S) 
ϕ  = Field(S, x->exp(4x))
D  = derivative(S)
I  = identity(S)

SU  = GaussLobatto{MU}(9)
SV  = GaussLobatto{MV}(9)
Sθ  = GaussLobatto{Mθ(9)

SUV  = SU ⦼ SV
SUVθ = SUV ⦼ Sθ
xUV  = Field(SUV)
xUVθ = Field(SUVθ)
Γ    = Field(SUV, (x,y)->x+y)

DI  = D ⦼ I
DII = ⦼(D, I, I)
DII = D ⊗ I ⊗ I
DII = (D ⊗ I) ⊗ I

#DEBUG
SUV = derivative(SUV)
