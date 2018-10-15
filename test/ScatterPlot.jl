using ScalarWave

struct U end
struct V end
struct UV end

#--------------------------------------------------------------------
# Define boundary and the product space
#--------------------------------------------------------------------
P1, P2 = 20, 20
M   = 1.0
SUV = ProductSpace{GaussLobatto{U,P1}, GaussLobatto{V,P2}}

#--------------------------------------------------------------------
# Define coordinates and their associated derivatives
#--------------------------------------------------------------------
𝕌 = Field(SUV, (U,V)->U)
𝕍 = Field(SUV, (U,V)->V)

𝑼 = ((-7M) + (-3M))/2 + (((-7M) - (-3M))/2)*𝕌  
𝑽 = (  7M  +   3M )/2 +  ( (7M  -   3M) /2)*𝕍  

t = Field(SUV, (𝑼,𝑽)->find_t_of_UV(𝑼, 𝑽, M), 𝑼, 𝑽)
r = Field(SUV, (𝑼,𝑽)->find_r_of_UV(𝑼, 𝑽, M), 𝑼, 𝑽)

ψ = 𝑼^2 + 𝑽^2

using PyPlot
fig = figure("pyplot_surfaceplot",figsize=(10,10))
ax  = fig[:add_subplot](2,1,1, projection = "3d")
ax[:plot_surface](t, r, ψ, rstride=2,edgecolors="k", cstride=2, cmap=ColorMap("gray"), alpha=0.8, linewidth=0.25)
xlabel("t")
ylabel("r")
title("Surface Plot")
show()
