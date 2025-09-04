using GLMakie, BeamletOptics

const BMO = BeamletOptics

GLMakie.activate!(; ssao=true)

include(joinpath(@__DIR__, "..", "render_utils.jl"))

mm = 1e-3
n = 1.5
cbs = RectangularPlateBeamsplitter(36mm, 25mm, 5mm, _->n)
cbs_mount = MeshDummy(joinpath(@__DIR__, "PBS Mount.stl"))

cbs_assembly = ObjectGroup([cbs, cbs_mount])

cmp = RectangularCompensatorPlate(36mm, 25mm, 5mm, _->n)
cmp_mount = MeshDummy(joinpath(@__DIR__, "PBS Mount.stl"))

cmp_assembly = ObjectGroup([cmp, cmp_mount])

zrotate3d!(cbs_assembly, deg2rad(45))
zrotate3d!(cmp_assembly, deg2rad(-45))

translate3d!(cmp_assembly, [0, 100mm, 0])

system = System([cbs_assembly, cmp_assembly])

beam = GaussianBeamlet([0,-100mm,0], [0, 11, 0], 1e-6, 5e-4)

solve_system!(system, beam)

##
pbs_view = [
    0.0940287   0.995569   -3.19189e-16  -0.048146
    -0.498165    0.0470503   0.865805      0.0176688
     0.861969   -0.0814105   0.500382     -0.155093
     0.0         0.0         0.0           1.0
]

fig = Figure(size=(600, 400)) 
display(fig)
ax = LScene(fig[1,1])
hide_axis(ax)

render!(ax, system)
render!(ax, beam)

set_view(ax, pbs_view)
save("pbs_showcase.png", fig; px_per_unit=8, update = false)