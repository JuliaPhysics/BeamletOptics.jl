using GLMakie, BeamletOptics
using LinearAlgebra

GLMakie.activate!(; ssao = true)

const BMO = BeamletOptics
const mm = 1e-3
const nm = 1e-9

##
rfl_val = 200mm
d_val = 60mm
angle_val = 90

oap = OffAxisParabolicMirror(rfl_val, d_val; angle = angle_val)
oap_sdf_obj = BMO.shape(oap)
f_parent_val = oap_sdf_obj.f
x_off_val = oap_sdf_obj.x_off

translate_to3d!(oap, [x_off_val, f_parent_val, 0])

beam = CollimatedSource(
    [x_off_val, -50mm, 0], [0, 1, 0], 25mm, 1550nm; num_rings = 4, num_rays = 80)

pd = Detector(60mm, true)
zrotate3d!(pd, deg2rad(angle_val))
translate_to3d!(pd, [0, f_parent_val, 0])

system = StaticSystem([oap, pd])
solve_system!(system, beam)

##
cview = [
  0.613294   0.789854  1.94289e-16  -0.156364
 -0.362033   0.281106  0.88877       0.0446795
  0.701999  -0.545077  0.458354     -0.222816
  0.0        0.0       0.0           1.0
]

fig = Figure(size = (600, 350))
display(fig)
ax = LScene(fig[1, 1], show_axis = false)

render!(ax, system)
render!(ax, beam)

set_view(ax, cview)

save("oap_mirror_showcase.png", fig; px_per_unit = 8, update = false)
