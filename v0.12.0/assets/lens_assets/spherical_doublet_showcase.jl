using GLMakie, BeamletOptics

GLMakie.activate!(; ssao=true)

const BMO = BeamletOptics
const mm = 1e-3
const cm = 10mm

include(joinpath(@__DIR__, "..", "render_utils.jl"))

##
λs = [488e-9, 707e-9, 1064e-9]

NLAK22 = DiscreteRefractiveIndex(λs, [1.6591, 1.6456, 1.6374])
NSF10 = DiscreteRefractiveIndex(λs, [1.7460, 1.7168, 1.7021])

AC254_150_AB = SphericalDoubletLens(87.9e-3, 105.6e-3, 1000, 6e-3, 3e-3, BeamletOptics.inch, NLAK22, NSF10)

system = System([AC254_150_AB])

##
cview = [
  0.602924   0.797798   0.0       -0.026697
 -0.108568   0.0820487  0.990697   0.000321103
  0.790377  -0.597315   0.136085  -0.0398802
  0.0        0.0        0.0        1.0
]

fig = Figure(size=(600,300))
ax = LScene(fig[1,1])

render!(ax, AC254_150_AB.front; alpha=0.3)
render!(ax, AC254_150_AB.back; alpha=0.3)

zs_1 = LinRange(-0.011, 0.011, 50)
# zs_2 = LinRange(-0.01, 0.01, 5)

for (i, z) in enumerate(zs_1)
    local beam = Beam([0, -0.02 , z], [0,1.,0], 488e-9)
    solve_system!(system, beam)
    render!(ax, beam, flen=0.15, color=RGBAf(0,0,1,0.7), show_pos=false)
end

for (i, z) in enumerate(zs_1)
    local beam = Beam([0, -0.02 , z], [0,1.,0], 707e-9)
    solve_system!(system, beam)
    render!(ax, beam, flen=0.15, color=RGBAf(1,0,0,0.5), show_pos=false)
end

display(fig)
hide_axis(ax)
set_view(ax, cview)
save("doublet_showcase.png", fig; px_per_unit=4, update = false)