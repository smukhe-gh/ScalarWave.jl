#--------------------------------------------------------------------
# Spacetime Discretization methods in Julia
# Soham 03-2018
#--------------------------------------------------------------------

using ScalarWave, Luxor

function backdrop(scene, framenumber)
    background("white")
end

struct U end
struct V end

function frame(scene, framenumber)
    field = Field(ProductSpace{GaussLobatto{U, 40}, GaussLobatto{V,40}},
                  (u,v) -> u*cos(rand()*pi*framenumber))
    drawpatch(field, "animation/animate-$(framenumber + 1000)")
end

function makeamovie()
    demo = Movie(800, 800, "test")
    animate(demo, [
        Scene(demo, backdrop, 0:300),
        Scene(demo, frame, 0:300, easingfunction=easeinoutcubic)
        ],
        creategif=true)
end

makeamovie()


