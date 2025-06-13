using GLMakie, BeamletOptics

GLMakie.activate!(; ssao=true)

const BMO = BeamletOptics

include(joinpath(@__DIR__, "..", "render_utils.jl"))

const mm = 1e-3

##
beam = Beam([0,0,0], [0,1,0], 1e-6)

m1 = RoundPlanoMirror(20mm, 5mm)
m2 = RoundPlanoMirror(20mm, 5mm)
m3 = RoundPlanoMirror(20mm, 5mm)

prism = RightAnglePrism(20mm, 20mm, n->1.2)

translate_to3d!(m1, [0, 100mm, 0])
zrotate3d!(m1, deg2rad(-30))

translate_to3d!(m2, [-70mm, 60mm, 0])
zrotate3d!(m2, deg2rad(180-75))

translate_to3d!(m3, [-35mm, 60mm, 0])
zrotate3d!(m3, deg2rad(180+45))

translate_to3d!(prism, [-28mm, 150mm, 0])
zrotate3d!(prism, deg2rad(-45))


system = System([m1, m2, m3, prism])

solve_system!(system, beam)

##
c_view = [
  0.0032104   0.999995    -5.34939e-5  -0.0948302
 -0.999993    0.00321049   0.00172712  -0.0287501
  0.00172728  4.79488e-5   0.999999    -6
  0.0         0.0          0.0          1.0
]

fig = Figure(size=(600,300))
display(fig)
ax = LScene(fig[1,1])
hide_axis(ax)
render!(ax, beam; flen=0, show_pos=true)
render!(ax, m1)
render!(ax, m2)
render!(ax, m3)
render!(ax, prism; color=RGBAf(1, 1, 1, .25))

arrow!(ax, position(first(BMO.rays(beam))), BMO.direction(first(BMO.rays(beam))); scale=2)


for _ray in BMO.rays(beam)
    dir = BMO.direction(_ray)
    pos = position(_ray)
    arrow!(ax, pos, dir; scale=2)
end

set_orthographic(ax)
set_view(ax, c_view)


save("beam_showcase.png", fig; px_per_unit=8, update = false)