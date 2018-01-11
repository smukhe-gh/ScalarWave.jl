#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 01-2018
#--------------------------------------------------------------------

# XXX: Also, occasionally test return types to check if your code 
# is type-stable.

using Base.Test
include("/Users/soham/Projects/spacetime-in-julia/src/SpecCalculus.jl") 

# test chebpoint
println("chebx: ", @test SpecCalculus.chebx(3, 10) ≈ 0.8090169943749475)

# test chebd
println("chebd: ", @test SpecCalculus.chebd(3, 1, 5) ≈ -0.723606797749979)

# test chebw
println("chebw: ", @test SpecCalculus.chebw(9, 12) ≈ 0.2258075258075258)