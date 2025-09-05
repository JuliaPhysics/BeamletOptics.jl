using GLMakie, BeamletOptics

GLMakie.activate!(; ssao=true)

const BMO = BeamletOptics
const mm = 1e-3

include(joinpath(@__DIR__, "..", "render_utils.jl"))

##
L3 = Lens(
    EvenAsphericalSurface(
        3.618e-3,               # r
        3.04e-3,                # d
        -44.874,                # conic
        [0,-0.14756*(1e3)^3, 0.035194*(1e3)^5, -0.0032262*(1e3)^7,
        0.0018592*(1e3)^9, 0.00036658*(1e3)^11, -0.00016039*(1e3)^13,
        -3.1846e-5*(1e3)^15]    # coeffs
    ),
    EvenAsphericalSurface(
        2.161e-3,               # r
        3.7e-3,                 # d
        -10.719,                # conic
        [0,-0.096568*(1e3)^3, 0.026771*(1e3)^5, -0.011261*(1e3)^7,
        0.0019879*(1e3)^9, 0.00015579*(1e3)^11, -0.00012433*(1e3)^13,
        1.5264e-5*(1e3)^15]     # coeffs
    ),
    0.7e-3,                     # center_thickness
    n -> 1.580200               # refractive index
)

xrotate3d!(L3, π/2)



##
cview = [
 -0.743072   0.669211  1.66533e-16   9.47028e-5
 -0.240303  -0.266825  0.933306      1.88076e-6
  0.624579   0.693513  0.359083     -0.00400377
  0.0        0.0       0.0           1.0
]
fig = Figure(size=(600, 300))
display(fig)
brightness = 1
ax = LScene(
    fig[1,1],
    scenekw = (
        lights = [
            DirectionalLight(RGBAf(brightness, brightness, brightness), Vec3f(-1, 0, 1))
            DirectionalLight(RGBAf(brightness, brightness, brightness), Vec3f(0, -1, -1))
        ], # hella weird syntax but whatever
    )
)
hide_axis(ax)

render!(ax, L3)
render!(ax, L3.shape.sdfs[3]; color=:blue, transparency=false)

set_view(ax, cview)
save("aspherical_lens_showcase.png", fig; px_per_unit=4, update = false)