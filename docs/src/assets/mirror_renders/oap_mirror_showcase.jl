using GLMakie, BeamletOptics
using LinearAlgebra

GLMakie.activate!(; ssao = true)

include(joinpath(@__DIR__, "..", "render_utils.jl"))

rfl_val = 200e-3
d_val = 60e-3
angle_val = 90

oap_mirror_obj = OffAxisParabolicMirror(rfl_val, d_val; angle = angle_val)
oap_sdf_obj = BeamletOptics.shape(oap_mirror_obj)
f_parent_val = oap_sdf_obj.f
x_off_val = oap_sdf_obj.x_off

translate_to3d!(oap_mirror_obj, [x_off_val, f_parent_val, 0])

beam_obj = CollimatedSource(
    [x_off_val, -50e-3, 0], [0, 1, 0], 25e-3, 1550e-9; num_rings = 4, num_rays = 80)

pd_obj = Detector(60e-3, true)
zrotate3d!(pd_obj, deg2rad(angle_val))
translate_to3d!(pd_obj, [0, f_parent_val, 0])

system_obj = BeamletOptics.StaticSystem([oap_mirror_obj, pd_obj])
solve_system!(system_obj, beam_obj)

fig_oap_sc = Figure(size = (600, 350))
ax_oap_sc = LScene(fig_oap_sc[1, 1], show_axis = false)

render!(ax_oap_sc, system_obj)
render!(ax_oap_sc, beam_obj)

cam3d!(ax_oap_sc.scene; eyeposition = Point3f(0.35, -0.15, 0.2),
    lookat = Point3f(0.1, 0.08, 0), upvector = Vec3f(0, 0, 1))

save("oap_mirror_showcase.png", fig_oap_sc; px_per_unit = 8, update = false)
