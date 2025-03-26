using CairoMakie, BeamletOptics

mm = 1e-3
n = 1.5
cbs = RectangularPlateBeamsplitter(36mm, 25mm, 5mm, _->n)
cbs_mount = MeshDummy(joinpath(file_dir, "PBS Mount.stl"))

cbs_assembly = ObjectGroup([cbs, cbs_mount])

cmp = RectangularCompensatorPlate(36mm, 25mm, 5mm, _->n)
cmp_mount = MeshDummy(joinpath(file_dir, "PBS Mount.stl"))

cmp_assembly = ObjectGroup([cmp, cmp_mount])

zrotate3d!(cbs_assembly, deg2rad(-45))
zrotate3d!(cmp_assembly, deg2rad(45))

translate3d!(cmp_assembly, [0, 100mm, 0])

system = System([cbs_assembly, cmp_assembly])

beam = GaussianBeamlet([0,-40mm,0], [0, 11, 0], 1e-6, 5e-4)

solve_system!(system, beam)

fig = Figure(size=(600, 320))
ax = Axis3(fig[1,1]; aspect=:data,  azimuth=0, elevation=pi/2)

render_object!(ax, cbs_mount)
render_object!(ax, cbs)
render_object!(ax, cmp_mount)
render_object!(ax, cmp)
render_beam!(ax, beam, flen=60mm, color=:red)

hidedecorations!(ax)
hidespines!(ax)

fp = ax.scene.plots[14]
ct = ax.scene.plots[15]
cm = ax.scene.plots[17]

fp.color[] = :white
ct.color[] = :magenta
cm.color[] = :white

fig