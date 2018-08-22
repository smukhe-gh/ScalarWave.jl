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

S  = GaussLobatto(9)
x  = Field(S) 
ϕ  = Field(S, x->exp(4x))
D  = derivative(S)

SU  = GaussLobatto(9)
SV  = GaussLobatto(9)
Sθ  = GaussLobatto(9)

SUV  = SU ⦼ SV
SUVθ = SUV ⦼ Sθ
xUV  = Field(SUV)
xUVθ = Field(SUVθ)
Γ    = Field(SUV, (x,y)->x+y)
