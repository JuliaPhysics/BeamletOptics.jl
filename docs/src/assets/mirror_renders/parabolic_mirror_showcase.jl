using GLMakie, BeamletOptics

GLMakie.activate!(; ssao = true)

const BMO = BeamletOptics
const mm = 1e-3
const nm = 1e-9

##
f_val = 100mm

pm = ParabolicMirror(f_val, 2BeamletOptics.inch)

beam = CollimatedSource([0, -200mm, 0], [0, 1, 0], 45mm, 1550nm; num_rings = 4, num_rays = 200)

system = StaticSystem([pm])
solve_system!(system, beam)

##
cview = [
  0.780432   0.625241  8.32667e-17   0.0589337
 -0.247203   0.308562  0.918521      0.0245199
  0.574297  -0.716843  0.395373     -0.133852
  0.0        0.0       0.0           1.0
]

fig = Figure(size = (600, 350))
display(fig)
ax = LScene(fig[1, 1], show_axis = false)

render!(ax, pm; transparency=true, alpha=0.5)
render!(ax, beam; flen = 2f_val, color=RGBAf(0,0,1,0.5), show_pos=false)

set_view(ax, cview)

save("parabolic_mirror_showcase.png", fig; px_per_unit = 8, update = false)
