#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 04-2018
# Test the derivative operators with ApproxFun and 
# ScalarWave
#--------------------------------------------------------------------

using ApproxFun, LinearAlgebra, HDF5

#--------------------------------------------------------------------
# Compute the action of the operator using ApproxFun
#--------------------------------------------------------------------
dU = -4 .. -3; 
dV =  3 ..  4;
d  = dU × dV
DU = ApproxFun.Derivative(d,[1,0]); 
DV = ApproxFun.Derivative(d,[0,1])

r  = Fun((U,V) -> find_r_of_UV(U,V,1.0), d)
invr = Fun((U,V) -> 1/find_r_of_UV(U,V,1.0), d)
UC = Fun((U,V) -> U, d)
VC = Fun((U,V) -> V, d)

DrDU = DU*r
DrDV = DV*r
DrDVDU = DU*DV*r

L  = DU*DV + ((DV*r)*invr)*DU + ((DU*r)*invr)*DV
ϕ  = (r^2*UC^2*VC)
Lϕ = L*ϕ

#--------------------------------------------------------------------
# Compute the action of the operator using ScalarWave 
#--------------------------------------------------------------------

M = 1.0
PV, PU = 29, 29
Umax, Umin = -3.0, -4.0
Vmin, Vmax =  3.0,  4.0
SUV = ScalarWave.ProductSpace{GaussLobatto(V,PV, Vmax, Vmin),
                   GaussLobatto(U,PU, Umax, Umin)}

𝔻𝕍, 𝔻𝕌 = derivative(SUV) 
𝕌  = Field(SUV, (U,V)->U)
𝕍  = Field(SUV, (U,V)->V)
𝕣  = Field(SUV, (U,V)->find_r_of_UV(U, V, M), 𝕌, 𝕍)
D𝕣D𝕌 = 𝔻𝕌*𝕣 
D𝕣D𝕍 = 𝔻𝕍*𝕣 
D𝕣D𝕍D𝕌 = 𝔻𝕌*𝔻𝕍*𝕣 

𝕃  = 𝔻𝕌*𝔻𝕍 + ((𝔻𝕌*𝕣)/𝕣)*𝔻𝕍 +((𝔻𝕍*𝕣)/𝕣)*𝔻𝕌
ψ  = (𝕣^2)*(𝕌^2)*𝕍
𝕃ψ = 𝕃*ψ

#--------------------------------------------------------------------
# Test components 
#--------------------------------------------------------------------

# test computation of r (𝕣)
r_array = zeros(30, 30)
DrDU_array = zeros(30, 30)
DrDV_array = zeros(30, 30)
DrDVDU_array = zeros(30, 30)

for _u in 1:30, _v in 1:30
    r_array[_u, _v] = r(𝕌.value[_u, _v], 𝕍.value[_u, _v]) 
end

for _u in 1:30, _v in 1:30
    DrDU_array[_u, _v] = DrDU(𝕌.value[_u, _v], 𝕍.value[_u, _v]) 
end

for _u in 1:30, _v in 1:30
    DrDV_array[_u, _v] = DrDV(𝕌.value[_u, _v], 𝕍.value[_u, _v]) 
end

for _u in 1:30, _v in 1:30
    DrDVDU_array[_u, _v] = DrDVDU(𝕌.value[_u, _v], 𝕍.value[_u, _v]) 
end

drawpatch(𝕣, "../output/scalar-wave-r")
drawpatch(Field(SUV, r_array), "../output/approxfun-r")

drawpatch(D𝕣D𝕌, "../output/scalarwave-drdu")
drawpatch(Field(SUV, DrDU_array), "../output/approxfun-drdu")

@test r_array ≈ 𝕣.value
@test DrDU_array ≈ D𝕣D𝕌.value
@test DrDV_array ≈ D𝕣D𝕍.value
@test DrDVDU_array ≈ D𝕣D𝕍D𝕌.value

