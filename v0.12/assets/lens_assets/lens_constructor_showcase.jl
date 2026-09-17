using GLMakie, BeamletOptics

GLMakie.activate!(; ssao=true)

const BMO = BeamletOptics
const mm = 1e-3
const inch = BMO.inch

include(joinpath(@__DIR__, "..", "render_utils.jl"))

##
r1 = 34.9e-3
r2 = -34.9e-3
l = 6.8e-3
LB1811 = Lens(
    SphericalSurface(r1, inch),
    SphericalSurface(r2, inch),
    l, 
    n -> 1.5
)

translate3d!(LB1811, [0, 10mm, 0])
system = System([LB1811])
b1 = Beam([0,-20mm,0], [0,1,0], 1e-6)
b2 = Beam([0,-20mm,5mm], [0,1,0], 1e-6)
b3 = Beam([0,-20mm,-5mm], [0,1,0], 1e-6)

solve_system!.(Ref(system), [b1, b2, b3])

##
fig = Figure(size=(600, 300))
display(fig)
ax = LScene(fig[1,1])
hide_axis(ax)
set_orthographic(ax)

render!(ax, LB1811; transparency=true, color=RGBf(0.678, 0.847, 0.902), alpha=0.25)
render!(ax, b1; show_pos=true, flen=33mm)
render!(ax, b2; show_pos=true, flen=33mm)
render!(ax, b3; show_pos=true, flen=33mm)

arrow!(ax, [0,0,0], [1,0,0]; color=:red)
arrow!(ax, [0,0,0], [0,1,0]; color=:green)
arrow!(ax, [0,0,0], [0,0,1]; color=:yellow)
text!(Point3(5mm, -2mm, 1mm), text="x")
text!(Point3(0, 5mm, 1mm), text="y")
text!(Point3(0, 0mm, 6.5mm), text="z")

cview = [
  0.010597    0.999943    -0.00118194  -0.0203364
 -0.0100706   0.00128867   0.999948    -0.000691592
  0.999893   -0.0105845    0.0100837   -1.63007
  0.0         0.0          0.0          1.0
]

set_view(ax, cview)
save("lens_constructor.png", fig; px_per_unit=4, update = false)