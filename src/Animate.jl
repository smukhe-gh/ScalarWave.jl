#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 08-2018
# Distorted Minkowski
#--------------------------------------------------------------------

addprocs(2)
@show nprocs() 

@everywhere using ScalarWave, Luxor

@everywhere struct U end
@everywhere struct V end
@everywhere struct UV end


@everywhere function frame(framenumber)
    P1, P2 = 60, 60
    SUV = ProductSpace{GaussLobatto{U,P1}, GaussLobatto{V,P2}}
    𝔹   = boundary(Null, SUV)
    
    u = Field(SUV, (u,v)->u)
    v = Field(SUV, (u,v)->v)
    Ω = Field(SUV, (u,v)->(framenumber*pi/120)*cospi(u/2)*cospi(v/2))
    
    𝒖 =  u*cos(Ω) + v*sin(Ω)
    𝒗 = -u*sin(Ω) + v*cos(Ω)
    𝔻𝒗, 𝔻𝒖 = derivativetransform(SUV, 𝒖, 𝒗)
    
    ρ = 0 
    𝕤 = exp(-((𝒖^2)/0.1)) 
    𝕓 = 𝔹*𝕤
    
    guu = Field(SUV, (u,v)-> 0)
    guv = Field(SUV, (u,v)->-2)
    gvv = Field(SUV, (u,v)-> 0)
    
    (𝕘𝒖𝒖, 𝕘𝒖𝒗, 𝕘𝒗𝒗) = inversemetrictransform(guu, guv, gvv, 𝒖, 𝒗) 
    invsqrtdet𝕘     = 1/sqrt(abs(inversemetricdet(𝕘𝒖𝒖, 𝕘𝒖𝒗, 𝕘𝒗𝒗))) 
    
    𝕘 = [𝕘𝒖𝒖 𝕘𝒖𝒗; 𝕘𝒖𝒗 𝕘𝒗𝒗]
    𝔻 = [𝔻𝒖, 𝔻𝒗]
    𝕃 = 𝕘𝒖𝒗*𝔻𝒖*𝔻𝒗 + 𝕘𝒖𝒗*𝔻𝒗*𝔻𝒖
    𝕨 = solve(𝕃 + 𝔹, ρ + 𝕓) 

    drawpatch(𝕨, "animate/minkowski-distorted-$(framenumber+1000)")
end


function makeamovie()
    @parallel vcat for framenumber in 1:30
        @show framenumber
        frame(framenumber)
    end
end

makeamovie()

