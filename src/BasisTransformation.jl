#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 09-2018
# Basis transformation functions FFT and PMMT
#--------------------------------------------------------------------

using ScalarWave, FFTW

order = 8

Λ = zeros(order+1, order+1)

for index in CartesianRange(size(Λ))
    n, m = (index.I .- 1)
    if m == 0 
        Λ[index] = 1/order
    elseif m == order
        Λ[index] = (1/order)*cospi(n)
    else
        Λ[index] = (2/order)*cospi((n*m)/order)
    end
end


𝕩 = chebgrid(order)
𝕔 = randn(order+1)  # 1D modal basis 
ℕ = [sum(𝕔[m+1]*cheb(m, x)for m in 0:order) for x in chebgrid(order)] 

#@show 𝕔 
@show (1/order)*(FFTW.r2r(ℕ, FFTW.REDFT00))./𝕔

𝕔 = randn(order+1, order+1)  # 2D modal basis 
ℕ = [sum(𝕔[m+1,n+1]*cheb(m, x)*cheb(n,y) for m in 0:order, n in 0:order) for x in chebgrid(order), y in chebgrid(order)] 

#@show 𝕔 
@show (1/order^2)*(FFTW.r2r(ℕ, FFTW.REDFT00))./𝕔

𝕔 = randn(order+1, order+1, order+1)  # 2D modal basis 
ℕ = [sum(𝕔[m+1,n+1,p+1]*cheb(m, x)*cheb(n,y)*cheb(p,z) 
         for m in 0:order, n in 0:order, p in 0:order) 
     for x in chebgrid(order), y in chebgrid(order), z in chebgrid(order)] 

#@show 𝕔 
@show (1/order^3)*(FFTW.r2r(ℕ, FFTW.REDFT00))./𝕔
