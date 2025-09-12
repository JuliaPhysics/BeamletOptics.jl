using GLMakie, BeamletOptics

GLMakie.activate!(; ssao=true)

const BMO = BeamletOptics

include(joinpath(@__DIR__, "..", "render_utils.jl"))

function spawn_mirror_mount()
    m1 = RoundPlanoMirror(1BeamletOptics.inch, 6e-3)
    holder = MeshDummy(joinpath(@__DIR__, "Mirror_Post.stl"))
    translate_to3d!(holder, [0,0,-(5.68e-2)])
    BeamletOptics.set_new_origin3d!(holder)
    mirror_mount = ObjectGroup([holder, m1])
    return mirror_mount
end

m1 = spawn_mirror_mount()
m2 = spawn_mirror_mount()
m3 = spawn_mirror_mount()
m4 = spawn_mirror_mount()

zrotate3d!(m1, deg2rad(45))
zrotate3d!(m2, deg2rad(-135))
zrotate3d!(m3, deg2rad(-45))
zrotate3d!(m4, deg2rad(135))

dy = 0.3

translate_to3d!(m1, [0,    0,      0])
translate_to3d!(m2, [0.1,  0,      0])
translate_to3d!(m3, [0.1,  dy,     0])
translate_to3d!(m4, [0,    dy,     0])

system = StaticSystem([m1, m2, m3, m4])

beam = GaussianBeamlet([0, -0.2, 0], [0, 1, 0], 532e-9, 1e-3)

solve_system!(system, beam, r_max=100)

## mirror render
mirror_camera = [
    0.723093   0.69075   2.498e-16  -0.075801
    -0.281368   0.294543  0.913278    0.0368553
    0.630847  -0.660385  0.407338   -0.254999
    0.0        0.0       0.0         1.0
]

fig = Figure(size=(600, 400)) 
display(fig)
ax = LScene(fig[1,1])
hide_axis(ax)

render!(ax, system)
render!(ax, beam; color=RGBf(0,1,0), flen=.4)

set_view(ax, mirror_camera)
save("plano_mirror_showcase.png", fig; px_per_unit=8, update = false)
