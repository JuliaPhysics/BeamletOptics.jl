using GLMakie, BeamletOptics

GLMakie.activate!(; ssao = true)

f_val = 100e-3

parabolic_mirror_obj = ParabolicMirror(f_val, 2BeamletOptics.inch)

beam_obj = CollimatedSource(
    [0, -150e-3, 0], [0, 1, 0], 20e-3, 1550e-9; num_rings = 4, num_rays = 80)

system_obj = StaticSystem([parabolic_mirror_obj])
solve_system!(system_obj, beam_obj)

fig_pm_sc = Figure(size = (600, 350))
ax_pm_sc = LScene(fig_pm_sc[1, 1], show_axis = false)

render!(ax_pm_sc, system_obj)
render!(ax_pm_sc, beam_obj; flen = 2f_val)

cam3d!(ax_pm_sc.scene; eyeposition = Point3f(0.05, -0.09, 0.2),
    lookat = Point3f(0, -0.09, 0), upvector = Vec3f(1, 0, 0))

save("parabolic_mirror_showcase.png", fig_pm_sc; px_per_unit = 8, update = false)
