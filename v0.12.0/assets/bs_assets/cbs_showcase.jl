using GLMakie, BeamletOptics

GLMakie.activate!(; ssao=true)

const BMO = BeamletOptics
const mm = 1e-3

include(joinpath(@__DIR__, "..", "render_utils.jl"))

##
n = 1.5

cbs1 = CubeBeamsplitter(BeamletOptics.inch, _->n)
cbs1_mount = MeshDummy(joinpath(@__DIR__, "CBS Mount.stl"))

cbs2 = CubeBeamsplitter(BeamletOptics.inch, _->n)
cbs2_mount = MeshDummy(joinpath(@__DIR__, "CBS Mount.stl"))
zrotate3d!(cbs2_mount, π)

cbs1_assembly = ObjectGroup([cbs1, cbs1_mount])
cbs2_assembly = ObjectGroup([cbs2, cbs2_mount])

m1 = RightAnglePrismMirror(BeamletOptics.inch, BeamletOptics.inch)
m1_mount = MeshDummy(joinpath(@__DIR__, "CBS Mount.stl"))
zrotate3d!(m1, π)

m2 = RightAnglePrismMirror(BeamletOptics.inch, BeamletOptics.inch)
m2_mount = MeshDummy(joinpath(@__DIR__, "CBS Mount.stl"))
zrotate3d!(m2, π)

m1_assembly = ObjectGroup([m1, m1_mount])
m2_assembly = ObjectGroup([m2, m2_mount])

translate_to3d!(cbs2_assembly, [-100mm, 150mm, 0])

translate_to3d!(m1_assembly, [-100mm, 0, 0])
zrotate3d!(m1_assembly, deg2rad(-45))

translate_to3d!(m2_assembly, [0, 150mm, 0])
zrotate3d!(m2_assembly, deg2rad(135))

system = System([cbs1_assembly, cbs2_assembly, m1_assembly, m2_assembly])

beam = GaussianBeamlet([0,-200mm,0], [0, 1, 0], 1e-6, 1e-3)

solve_system!(system, beam)

## splitter render
mzi_view = [
    -0.691546   0.722333  -5.27356e-16  -0.0997083
    -0.540886  -0.517833   0.662791      0.0538701
    0.478755   0.45835    0.748805     -0.3149
    0.0        0.0        0.0           1.0
]

fig = Figure(size=(600, 400)) 
display(fig)
ax = LScene(fig[1,1])
hide_axis(ax)

render!(ax, system)
render!(ax, beam)

set_view(ax, mzi_view)
save("cbs_showcase.png", fig; px_per_unit=8, update = false)