#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 09-2018
# Basis transformation functions FFT and PMMT
#--------------------------------------------------------------------

using ScalarWave, FFTW

order = 7
𝕩 = chebgrid(order)
𝕔 = randn(order+1)  # 1D modal basis 
𝕟 = [sum(𝕔[m+1]*cheb(m, x)for m in 0:order) for x in chebgrid(order)] 

dct𝕔  = (1/order)*(FFTW.r2r(𝕟, FFTW.REDFT00))
idct𝕔 = (FFTW.r2r(order*dct𝕔, FFTW.REDFT00))/(2*order)
@show 𝕟
@show 𝕔
@show dct𝕔./𝕔
@show idct𝕔./𝕟
